"""
10_hiv_history.py — HIV disease-history factors and 2025 seropositivity (within PWH)

Outcome: anti-PT IgG >=10 EU/mL in 2025. Population: PWH 2025 (frozen n = 528).
Exposures come from the baseline record sheet of the 2021 study (SRC_2021_CSV),
linked to the 2025 cohort on participant ID + date of birth:
  - CD4 count at ART initiation (recorded as such; not a true nadir)
  - ART duration = ART start date to the 2025 sampling date
  - AIDS history (yes/no)
  - current CD4 count (2025, from the analysis master)

ART start dates are hand-entered and heterogeneous. Parsing rules:
  - YYYY/M/D ............ as given
  - 'Oct-09', 'YYYY/M/' . month known -> 15th of that month
  - 'YYYY', 'YYYY/' ..... year only -> 1 July (+/- 6 months)
  - 5-digit number ...... Excel serial date (origin 1899-12-30)
  - elite controller .... never started ART -> ART duration missing (counted separately)
  - '/' or blank ........ missing
The resulting error (at most +/- 0.5 years) does not affect the categories
(<10 / 10-19 / >=20 years) or the 5-year continuous term.

Added analyses (within PWH 2025):
  - CD4/CD8 ratio at the 2025 sampling (CD8 from SRC_CD8_2025_CSV, linked on ID + date of
    birth), <=0.5 vs >0.5, against 2025 seropositivity -> numbers/cd4cd8_2025.csv.
    Read together with the 2021 baseline ratio (06).
  - Priming-vaccine generation by birth year (<=1974 wP / 1975-81 transition / >=1982 aP;
    same boundaries as 09) -> numbers/pwh_priming_generation.csv. Contrast with PWoH (09).

Input : output/analysis_master.csv, SRC_2021_CSV, SRC_CD8_2025_CSV
Output: numbers/hiv_history.csv (univariable), numbers/hiv_history_model.csv (multivariable),
        numbers/cd4cd8_2025.csv, numbers/pwh_priming_generation.csv,
        output/hiv_history_cc.csv (per-participant input for verify_in_stata.do; private)
"""
import re
import sys

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

import config as C

CD4_ART_BINS = [-np.inf, 199, 349, np.inf]
CD4_ART_LABELS = ["<200", "200-349", ">=350"]
ART_YR_BINS = [-np.inf, 9.999, 19.999, np.inf]
ART_YR_LABELS = ["<10", "10-19", ">=20"]
CD4_NOW_BINS = [-np.inf, 349, 499, np.inf]
CD4_NOW_LABELS = ["<350", "350-499", ">=500"]


def parse_art_start(v) -> tuple:
    """Parse a heterogeneous ART-start entry. Returns (Timestamp or None, kind)."""
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return None, "missing"
    s = str(v).strip()
    if s in ("", "/", "nan"):
        return None, "missing"
    if "コントローラー" in s:          # "controller": never started ART
        return None, "controller"
    if re.fullmatch(r"\d{5}(\.0)?", s):
        return pd.Timestamp("1899-12-30") + pd.Timedelta(days=int(float(s))), "serial"
    m = re.fullmatch(r"(\d{4})[/-](\d{1,2})[/-](\d{1,2})", s)
    if m:
        return pd.Timestamp(int(m[1]), int(m[2]), int(m[3])), "date"
    m = re.fullmatch(r"(\d{4})[/-](\d{1,2})[/-]?", s)
    if m:
        return pd.Timestamp(int(m[1]), int(m[2]), 15), "month"
    m = re.fullmatch(r"([A-Za-z]{3})-(\d{2})", s)
    if m:
        mon = pd.to_datetime(m[1], format="%b").month
        return pd.Timestamp(2000 + int(m[2]), mon, 15), "month"
    m = re.fullmatch(r"(\d{4})(\.0)?[/?]?", s)
    if m:
        return pd.Timestamp(int(m[1]), 7, 1), "year"
    return None, "unparsed"


def to_cd4(v):
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return np.nan
    s = re.sub(r"[^0-9.]", "", str(v))   # strips qualifiers such as 'approx.'
    return float(s) if s else np.nan


def load_history() -> pd.DataFrame:
    d = pd.read_csv(C.SRC_2021_CSV).rename(columns={
        "生年月日": "dob_raw", "IMS no.": "id_raw", "ART開始日": "art_raw",
        "AIDS歴": "aids_raw", "ART開始時のCD4数": "cd4_art_raw",
    })
    d["dob"] = pd.to_datetime(d["dob_raw"], format="mixed")
    d["id"] = (d["id_raw"].astype(str).str.upper()
               .str.replace("-", "", regex=False).str.strip())
    parsed = d["art_raw"].map(parse_art_start)
    d["art_start"] = [p[0] for p in parsed]
    d["art_parse"] = [p[1] for p in parsed]
    d["aids"] = d["aids_raw"].map({"有": 1, "無": 0})   # yes / no
    d["cd4_art"] = d["cd4_art_raw"].map(to_cd4)
    d["art_naive"] = (d["art_parse"] == "controller").astype(int)
    key = d["id"] + "|" + d["dob"].dt.strftime("%Y%m%d")
    if key.duplicated().any():
        raise AssertionError("2021 record sheet: duplicate ID + date of birth")
    return d[["id", "dob", "art_start", "art_parse", "art_naive", "aids", "cd4_art"]]


def load_cd8_2025() -> pd.DataFrame:
    """CD8 count at the 2025 sampling (de-identified extract), keyed on ID + date of birth."""
    d = pd.read_csv(C.SRC_CD8_2025_CSV, parse_dates=["dob"])
    d = d.dropna(subset=["id"]).copy()
    d["id"] = d["id"].astype(str).str.upper().str.replace("-", "", regex=False).str.strip()
    d = d.drop_duplicates(subset=["id", "dob"])
    return d[["id", "dob", "cd8"]].rename(columns={"cd8": "cd8_2025"})


def gmc(s):
    x = np.log(s.astype(float)); n = len(x)
    se = x.std(ddof=1) / np.sqrt(n); t = stats.t.ppf(0.975, n - 1)
    return np.exp(x.mean()), np.exp(x.mean() - t * se), np.exp(x.mean() + t * se)


def univariate(df, var, labels, name):
    rows = []
    sub = df[df[var].notna()]
    groups = [sub[sub[var] == lab] for lab in labels]
    tab = np.array([[g["pt_pos"].sum(), len(g) - g["pt_pos"].sum()] for g in groups])
    if len(labels) == 2:
        p_pos = stats.fisher_exact(tab)[1]; test = "Fisher"
        p_titre = stats.mannwhitneyu(groups[0]["pt"], groups[1]["pt"],
                                     alternative="two-sided").pvalue; ttest = "MWU"
    else:
        p_pos = stats.chi2_contingency(tab, correction=False)[1]; test = "chi2"
        p_titre = stats.kruskal(*[g["pt"] for g in groups]).pvalue; ttest = "KW"
    for lab, g in zip(labels, groups):
        gm = gmc(g["pt"])
        rows.append({"variable": name, "level": lab, "N": len(g),
                     "seropositive_n": int(g["pt_pos"].sum()),
                     "seropositive_pct": round(100 * g["pt_pos"].mean(), 1),
                     "GMC": round(gm[0], 1), "GMC_lci": round(gm[1], 1), "GMC_uci": round(gm[2], 1),
                     "PT_median": g["pt"].median()})
    rows[0].update({"test_pos": test, "P_pos": round(float(p_pos), 4),
                    "test_titre": ttest, "P_titre": round(float(p_titre), 4),
                    "N_missing": int(len(df) - len(sub))})
    return rows, float(p_pos), float(p_titre)


def main():
    print("=" * 72)
    print("10_hiv_history — HIV disease history x 2025 seropositivity (within PWH)")
    print("=" * 72)

    m = pd.read_csv(C.MASTER, parse_dates=["dob", "draw_date"])
    p = m[m["hiv"] == 1].copy()
    C.check("n_pwh_2025", len(p), label="PWH 2025")

    h = load_history()
    print("  ART-start entry types: " + ", ".join(f"{k}={v}" for k, v in h["art_parse"].value_counts().items()))
    p = p.merge(h, on=["id", "dob"], how="left")
    assert len(p) == C.EXPECTED["n_pwh_2025"], "row count changed on merge"
    linked = p["aids"].notna() | p["cd4_art"].notna() | p["art_start"].notna()
    print(f"  PWH linked to the 2021 record sheet: {int(linked.sum())}/{len(p)}")

    # 2025 CD8 -> 2025 CD4/CD8 ratio
    c8 = load_cd8_2025()
    p = p.merge(c8, on=["id", "dob"], how="left")
    assert len(p) == C.EXPECTED["n_pwh_2025"], "row count changed on CD8 merge"
    p["cd48_ratio_2025"] = p["cd4"] / p["cd8_2025"]
    print(f"  2025 CD8 available: {int(p['cd8_2025'].notna().sum())}/{len(p)}")

    # priming-vaccine generation (same boundaries as 09)
    by = p["dob"].dt.year
    p["generation"] = np.where(by <= 1974, "wP", np.where(by >= 1982, "aP", "transition"))

    p["art_years"] = (p["draw_date"] - p["art_start"]).dt.days / 365.25
    p["cd4_art_cat"] = pd.cut(p["cd4_art"], CD4_ART_BINS, labels=CD4_ART_LABELS)
    p["art_cat"] = pd.cut(p["art_years"], ART_YR_BINS, labels=ART_YR_LABELS)
    p["cd4_now_cat"] = pd.cut(p["cd4"], CD4_NOW_BINS, labels=CD4_NOW_LABELS)
    p["aids_cat"] = p["aids"].map({0: "no", 1: "yes"})
    n_controller = int(p["art_naive"].fillna(0).sum())
    print(f"  ART-naive (controller): {n_controller} -> ART duration treated as missing")

    # per-participant file for the Stata recomputation (no identifiers; output/ is private)
    cols = ["pt_pos", "pt", "age", "male", "cd4", "cd4_art", "art_years", "aids",
            "cd4_now_cat", "cd4_art_cat", "art_cat", "cd8_2025", "cd48_ratio_2025", "generation"]
    p[cols].to_csv(C.DIR_OUT / "hiv_history_cc.csv", index=False)

    rows = []; pvals = {}
    for var, labels, name in [("cd4_now_cat", CD4_NOW_LABELS, "CD4 2025 (cells/uL)"),
                              ("cd4_art_cat", CD4_ART_LABELS, "CD4 at ART initiation (cells/uL)"),
                              ("art_cat", ART_YR_LABELS, "ART duration at 2025 sampling (years)"),
                              ("aids_cat", ["no", "yes"], "AIDS history")]:
        r, pp, pt = univariate(p, var, labels, name)
        rows += r; pvals[var] = (pp, pt)
        print(f"  {name}: " + " / ".join(f"{x['level']} {x['seropositive_pct']}% (n={x['N']})" for x in r)
              + f"  P={pp:.3f} (titre P={pt:.3f})")
    pd.DataFrame(rows).to_csv(C.DIR_NUM / "hiv_history.csv", index=False)

    # multivariable model, complete cases within PWH
    cc = p.dropna(subset=["age", "male", "cd4", "cd4_art", "art_years", "aids"]).copy()
    X = sm.add_constant(pd.DataFrame({
        "age_per10": cc["age"] / 10,
        "male": cc["male"].astype(float),
        "cd4_2025_per100": cc["cd4"] / 100,
        "cd4_art_per100": cc["cd4_art"] / 100,
        "art_per5y": cc["art_years"] / 5,
        "aids": cc["aids"].astype(float),
    }))
    fit = sm.Logit(cc["pt_pos"].astype(float), X).fit(disp=0)
    ci = fit.conf_int()
    mod = pd.DataFrame({
        "term": fit.params.index, "aOR": np.exp(fit.params.values).round(3),
        "lci": np.exp(ci[0].values).round(3), "uci": np.exp(ci[1].values).round(3),
        "P": fit.pvalues.values.round(4)})
    mod = mod[mod["term"] != "const"]
    mod.insert(0, "N", len(cc))
    mod.to_csv(C.DIR_NUM / "hiv_history_model.csv", index=False)
    print(f"\n  multivariable (n={len(cc)}):\n" + mod.to_string(index=False))

    # 2025 CD4/CD8 ratio (<=0.5 / >0.5) x 2025 seropositivity
    r = p.dropna(subset=["cd48_ratio_2025"]).copy()
    r["stratum"] = np.where(r["cd48_ratio_2025"] <= C.CD48_CUT, "low", "high")
    lo = r[r["stratum"] == "low"]; hi = r[r["stratum"] == "high"]
    tab = np.array([[lo["pt_pos"].sum(), len(lo) - lo["pt_pos"].sum()],
                    [hi["pt_pos"].sum(), len(hi) - hi["pt_pos"].sum()]])
    or25, p25 = stats.fisher_exact(tab)
    p25_mwu = stats.mannwhitneyu(lo["pt"], hi["pt"], alternative="two-sided").pvalue
    Xr = sm.add_constant(pd.DataFrame({"age_per10": r["age"] / 10, "male": r["male"].astype(float),
                                       "ratio": r["cd48_ratio_2025"]}))
    fr = sm.Logit(r["pt_pos"].astype(float), Xr).fit(disp=0)
    cir = np.exp(fr.conf_int().loc["ratio"])
    g_lo = gmc(lo["pt"]); g_hi = gmc(hi["pt"])
    out25 = pd.DataFrame([
        {"stratum": f"CD4/CD8 2025 <={C.CD48_CUT}", "N": len(lo), "seropositive_n": int(lo["pt_pos"].sum()),
         "seropositive_pct": round(100 * lo["pt_pos"].mean(), 1), "GMC": round(g_lo[0], 1),
         "GMC_lci": round(g_lo[1], 1), "GMC_uci": round(g_lo[2], 1), "PT_median": lo["pt"].median()},
        {"stratum": f"CD4/CD8 2025 >{C.CD48_CUT}", "N": len(hi), "seropositive_n": int(hi["pt_pos"].sum()),
         "seropositive_pct": round(100 * hi["pt_pos"].mean(), 1), "GMC": round(g_hi[0], 1),
         "GMC_lci": round(g_hi[1], 1), "GMC_uci": round(g_hi[2], 1), "PT_median": hi["pt"].median()},
    ])
    out25["OR_fisher"] = [round(float(or25), 3), ""]
    out25["P_fisher"] = [round(float(p25), 4), ""]
    out25["P_mwu"] = [round(float(p25_mwu), 4), ""]
    out25["aOR_per_unit_age_sex_adj"] = [round(float(np.exp(fr.params["ratio"])), 3), ""]
    out25["aOR_lci"] = [round(float(cir[0]), 3), ""]; out25["aOR_uci"] = [round(float(cir[1]), 3), ""]
    out25["aOR_P"] = [round(float(fr.pvalues["ratio"]), 4), ""]
    out25["ratio_median_2025"] = [round(float(r["cd48_ratio_2025"].median()), 2), ""]
    both = r.dropna(subset=["cd48_ratio_2021"])
    rho = stats.spearmanr(both["cd48_ratio_2021"], both["cd48_ratio_2025"]).statistic
    out25["spearman_2021_vs_2025"] = [round(float(rho), 3), ""]; out25["n_both_ratios"] = [len(both), ""]
    out25.to_csv(C.DIR_NUM / "cd4cd8_2025.csv", index=False)
    print(f"\n  2025 CD4/CD8 ratio: <=0.5 n={len(lo)} {100*lo['pt_pos'].mean():.1f}% / >0.5 n={len(hi)} "
          f"{100*hi['pt_pos'].mean():.1f}%  Fisher p={p25:.3f}, MWU p={p25_mwu:.3f}, "
          f"aOR/unit {np.exp(fr.params['ratio']):.2f} (p={fr.pvalues['ratio']:.3f}); "
          f"Spearman rho vs 2021 ratio={rho:.2f} (n={len(both)})")
    C.check("cd48_2025_n", len(r), label="2025 CD4/CD8 n")
    C.check("cd48_2025_low_n", len(lo), label="2025 ratio <=0.5 N")
    C.check("cd48_2025_high_n", len(hi), label="2025 ratio >0.5 N")
    C.check_pct("cd48_2025_low_pos_pct", 100 * lo["pt_pos"].mean(), label="2025 ratio <=0.5 seropositive %")
    C.check_pct("cd48_2025_high_pos_pct", 100 * hi["pt_pos"].mean(), label="2025 ratio >0.5 seropositive %")
    C.check("cd48_2025_fisher_p", float(p25), label="2025 ratio Fisher p")
    C.check("cd48_2025_mwu_p", float(p25_mwu), label="2025 ratio MWU p")

    # priming-vaccine generation within PWH
    rows_g = []
    for gname in ("wP", "transition", "aP"):
        g = p[p["generation"] == gname]
        gm = gmc(g["pt"])
        rows_g.append({"generation": gname, "N": len(g), "seropositive_n": int(g["pt_pos"].sum()),
                       "seropositive_pct": round(100 * g["pt_pos"].mean(), 1),
                       "GMC": round(gm[0], 1), "GMC_lci": round(gm[1], 1), "GMC_uci": round(gm[2], 1)})
    ga = p[p["generation"] == "aP"]; gw = p[p["generation"] == "wP"]
    tabg = np.array([[ga["pt_pos"].sum(), len(ga) - ga["pt_pos"].sum()],
                     [gw["pt_pos"].sum(), len(gw) - gw["pt_pos"].sum()]])
    org, pg = stats.fisher_exact(tabg)
    outg = pd.DataFrame(rows_g)
    outg["OR_aP_vs_wP"] = [round(float(org), 3), "", ""]
    outg["P_fisher_aP_vs_wP"] = [round(float(pg), 4), "", ""]
    outg.to_csv(C.DIR_NUM / "pwh_priming_generation.csv", index=False)
    print("\n  PWH by priming generation:\n" + outg.to_string(index=False))
    C.check("pwh_gen_ap_n", len(ga), label="PWH aP N")
    C.check("pwh_gen_wp_n", len(gw), label="PWH wP N")
    C.check("pwh_gen_trans_n", int((p["generation"] == "transition").sum()), label="PWH transition N")
    C.check_pct("pwh_gen_ap_pos_pct", 100 * ga["pt_pos"].mean(), label="PWH aP seropositive %")
    C.check_pct("pwh_gen_wp_pos_pct", 100 * gw["pt_pos"].mean(), label="PWH wP seropositive %")
    C.check("pwh_gen_fisher_p", float(pg), label="PWH aP vs wP Fisher p")

    C.check("hivhx_linked_n", int(linked.sum()), label="linked to 2021 sheet")
    C.check("hivhx_model_n", len(cc), label="HIV-history model n")
    C.check("hivhx_cd4art_p", pvals["cd4_art_cat"][0], label="CD4 at ART P")
    C.check("hivhx_artdur_p", pvals["art_cat"][0], label="ART duration P")
    C.check("hivhx_aids_p", pvals["aids_cat"][0], label="AIDS history P")
    C.check("hivhx_cd4now_p", pvals["cd4_now_cat"][0], label="CD4 2025 P")
    print("  -> numbers/hiv_history.csv, numbers/hiv_history_model.csv, numbers/cd4cd8_2025.csv, "
          "numbers/pwh_priming_generation.csv")

    return {"hivhx_linked_n": int(linked.sum()), "hivhx_model_n": len(cc),
            "cd48_2025": {"n": len(r), "low_pct": round(100 * lo["pt_pos"].mean(), 1),
                          "high_pct": round(100 * hi["pt_pos"].mean(), 1), "fisher_p": round(float(p25), 4),
                          "mwu_p": round(float(p25_mwu), 4), "aOR_per_unit": round(float(np.exp(fr.params["ratio"])), 3),
                          "aOR_p": round(float(fr.pvalues["ratio"]), 4), "spearman_vs_2021": round(float(rho), 3)},
            "pwh_generation": {"aP_pct": round(100 * ga["pt_pos"].mean(), 1), "wP_pct": round(100 * gw["pt_pos"].mean(), 1),
                               "fisher_p": round(float(pg), 4)},
            "hivhx_controller_n": n_controller,
            "hivhx_univariate_p": {k: round(v[0], 4) for k, v in pvals.items()},
            "hivhx_model": {t: {"aOR": a, "lci": l, "uci": u, "P": pp}
                            for t, a, l, u, pp in mod[["term", "aOR", "lci", "uci", "P"]].itertuples(index=False)}}


if __name__ == "__main__":
    main()
    sys.exit(0)
