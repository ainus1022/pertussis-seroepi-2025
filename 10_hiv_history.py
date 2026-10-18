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

Input : output/analysis_master.csv, SRC_2021_CSV
Output: numbers/hiv_history.csv (univariable), numbers/hiv_history_model.csv (multivariable),
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

    p["art_years"] = (p["draw_date"] - p["art_start"]).dt.days / 365.25
    p["cd4_art_cat"] = pd.cut(p["cd4_art"], CD4_ART_BINS, labels=CD4_ART_LABELS)
    p["art_cat"] = pd.cut(p["art_years"], ART_YR_BINS, labels=ART_YR_LABELS)
    p["cd4_now_cat"] = pd.cut(p["cd4"], CD4_NOW_BINS, labels=CD4_NOW_LABELS)
    p["aids_cat"] = p["aids"].map({0: "no", 1: "yes"})
    n_controller = int(p["art_naive"].fillna(0).sum())
    print(f"  ART-naive (controller): {n_controller} -> ART duration treated as missing")

    # per-participant file for the Stata recomputation (no identifiers; output/ is private)
    cols = ["pt_pos", "pt", "age", "male", "cd4", "cd4_art", "art_years", "aids",
            "cd4_now_cat", "cd4_art_cat", "art_cat"]
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

    C.check("hivhx_linked_n", int(linked.sum()), label="linked to 2021 sheet")
    C.check("hivhx_model_n", len(cc), label="HIV-history model n")
    C.check("hivhx_cd4art_p", pvals["cd4_art_cat"][0], label="CD4 at ART P")
    C.check("hivhx_artdur_p", pvals["art_cat"][0], label="ART duration P")
    C.check("hivhx_aids_p", pvals["aids_cat"][0], label="AIDS history P")
    C.check("hivhx_cd4now_p", pvals["cd4_now_cat"][0], label="CD4 2025 P")
    print("  -> numbers/hiv_history.csv, numbers/hiv_history_model.csv")

    return {"hivhx_linked_n": int(linked.sum()), "hivhx_model_n": len(cc),
            "hivhx_controller_n": n_controller,
            "hivhx_univariate_p": {k: round(v[0], 4) for k, v in pvals.items()},
            "hivhx_model": {t: {"aOR": a, "lci": l, "uci": u, "P": pp}
                            for t, a, l, u, pp in mod[["term", "aOR", "lci", "uci", "P"]].itertuples(index=False)}}


if __name__ == "__main__":
    main()
    sys.exit(0)
