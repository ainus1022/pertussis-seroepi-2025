"""
02_table1_characteristics.py - Table 1 and seropositivity by age group (PWoH vs PWH).

Input : output/analysis_master.csv, output/frozen_biobank.csv   (private)
Output: numbers/table1.csv, numbers/table1_agestrat.csv,
        numbers/table1_agestrat_setB.csv, figures/Fig1_*

PWoH is shown as two sets (config section 10):
  A = hospital outpatients (301) + serum bank (149) = 450   main comparison
  B = hospital outpatients only (301)                       has blood counts and reason for visit
Variables that do not exist for the serum bank are shown as "-" for set A.

The primary outcome is seropositivity (anti-PT IgG >= 10 EU/mL); geometric mean
concentrations are reported as supporting information.
"""
import sys
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import config as C
import lib

CALLER = "02_table1_characteristics"

AGE_BINS = [0, 40, 50, 60, 70, 200]
AGE_LABELS = ["<40", "40-49", "50-59", "60-69", "70+"]

NA = "-"   # structurally unavailable for set A


def med_iqr(s: pd.Series, dec: int = 0) -> str:
    s = s.dropna()
    if s.empty:
        return NA
    return f"{s.median():.{dec}f} [{s.quantile(.25):.{dec}f}-{s.quantile(.75):.{dec}f}]"


def gmc(s: pd.Series):
    """Geometric mean and 95% CI (t distribution on the log scale)."""
    x = np.log(s.dropna().astype(float))
    n = len(x)
    se = x.std(ddof=1) / np.sqrt(n)
    t = stats.t.ppf(0.975, n - 1)
    return np.exp(x.mean()), np.exp(x.mean() - t * se), np.exp(x.mean() + t * se)


def p_mwu(a: pd.Series, b: pd.Series):
    a, b = a.dropna(), b.dropna()
    if a.empty or b.empty:
        return None
    return float(stats.mannwhitneyu(a, b).pvalue)


def p_chi2(k1, n1, k0, n0):
    tab = np.array([[k1, n1 - k1], [k0, n0 - k0]])
    return float(stats.chi2_contingency(tab)[1])


def fmt_p(p):
    if p is None:
        return NA
    return "<0.001" if p < 0.001 else f"{p:.3g}"


def age_strat(pwoh: pd.DataFrame, pwh: pd.DataFrame) -> pd.DataFrame:
    """Seropositivity by age group, PWoH vs PWH, with Fisher's exact test per stratum."""
    def binned(d):
        d = d.copy()
        d["age_grp"] = pd.cut(d["age"], bins=AGE_BINS, labels=AGE_LABELS, right=False)
        return d

    a, p = binned(pwoh), binned(pwh)
    rows = []
    for g in AGE_LABELS:
        sa, sp = a[a["age_grp"] == g], p[p["age_grp"] == g]
        if len(sa) == 0 or len(sp) == 0:
            continue
        ka, kp = int(sa["pt_pos"].sum()), int(sp["pt_pos"].sum())
        tab = np.array([[ka, len(sa) - ka], [kp, len(sp) - kp]])
        rows.append({
            "age_grp": g,
            "PWoH_N": len(sa), "PWoH_pos": ka, "PWoH_pct": round(100 * ka / len(sa), 1),
            "PWH_N": len(sp), "PWH_pos": kp, "PWH_pct": round(100 * kp / len(sp), 1),
            "P_fisher": round(float(stats.fisher_exact(tab)[1]), 4),
        })
    return pd.DataFrame(rows)


def main():
    print("=" * 72)
    print("02_table1_characteristics - Table 1 (PWoH sets A and B)")
    print("=" * 72)

    m = pd.read_csv(C.MASTER)
    pwh = m[m["hiv"] == 1]
    A = lib.load_pwoh_set(CALLER, "A")     # 450
    B = lib.load_pwoh_set(CALLER, "B")     # 301

    C.check("pwoh_a_n", len(A), label="PWoH-A n")
    C.check("n_pwoh_2025", len(B), label="PWoH-B n")
    print(f"  PWH {len(pwh)} / PWoH-A {len(A)} (outpatients {int((A['cohort'] == 'main').sum())} "
          f"+ serum bank {int((A['cohort'] == 'biobank').sum())}) / PWoH-B {len(B)}")

    cPWH = f"PWH (n={len(pwh)})"
    cA = f"PWoH-A (n={len(A)})"
    cB = f"PWoH-B (n={len(B)})"
    rows = []

    def add(label, v_pwh, v_a, v_b, p_a="", p_b=""):
        rows.append({"Variable": label, cPWH: v_pwh, cA: v_a, cB: v_b,
                     "P (PWH vs A)": p_a, "P (PWH vs B)": p_b})

    # --- Age ---
    add("Age, years - median [IQR]",
        med_iqr(pwh["age"]), med_iqr(A["age"]), med_iqr(B["age"]),
        fmt_p(p_mwu(pwh["age"], A["age"])), fmt_p(p_mwu(pwh["age"], B["age"])))
    C.check("pwh_age_median", pwh["age"].median(), label="PWH median age")
    C.check("pwoh_age_median", B["age"].median(), label="PWoH-B median age")
    C.check("pwoh_a_age_median", A["age"].median(), label="PWoH-A median age")
    C.check("table1_a_age_p", p_mwu(pwh["age"], A["age"]), label="PWH vs A, age p")

    # --- Sex ---
    def sex_cell(d):
        return f"{int(d['male'].sum())} ({100 * d['male'].mean():.1f})"
    add("Male - n (%)", sex_cell(pwh), sex_cell(A), sex_cell(B),
        fmt_p(p_chi2(pwh["male"].sum(), len(pwh), A["male"].sum(), len(A))),
        fmt_p(p_chi2(pwh["male"].sum(), len(pwh), B["male"].sum(), len(B))))
    C.check_pct("pwh_male_pct", 100 * pwh["male"].mean(), label="PWH male %")
    C.check_pct("pwoh_a_male_pct", 100 * A["male"].mean(), label="PWoH-A male %")

    # --- Blood counts (not available for set A) ---
    add("WBC, /uL - median [IQR]", med_iqr(pwh["wbc"]), NA, med_iqr(B["wbc"]),
        NA, fmt_p(p_mwu(pwh["wbc"], B["wbc"])))
    add("Lymphocytes, /uL - median [IQR]", med_iqr(pwh["lym"]), NA, med_iqr(B["lym"]),
        NA, fmt_p(p_mwu(pwh["lym"], B["lym"])))
    C.check("pwh_wbc_median", pwh["wbc"].median(), label="PWH median WBC")
    C.check("pwh_lym_median", pwh["lym"].median(), label="PWH median lymphocytes")
    C.check("pwoh_lym_median", B["lym"].median(), label="PWoH-B median lymphocytes")

    # --- HIV-related (PWH only) ---
    add("CD4, /uL - median [IQR]", med_iqr(pwh["cd4"]), NA, NA)
    C.check("pwh_cd4_median", pwh["cd4"].median(), label="PWH median CD4")
    vl = pwh["vl"].dropna()
    n_supp = int((vl < 50).sum())
    add("HIV-RNA <50 copies/mL - n/N (%)",
        f"{n_supp}/{len(vl)} ({100 * n_supp / len(vl):.1f})", NA, NA)
    C.check("pwh_vl_supp_n", n_supp, label="HIV-RNA <50 n")
    C.check("pwh_vl_supp_denom", len(vl), label="HIV-RNA measured n")

    # --- Reason for visit (set B only; not applicable to serum-bank volunteers) ---
    cat = B["dx_cat"].value_counts()
    for key, lab in [("liver", "  Liver disease"), ("metabolic", "  Metabolic disease"),
                     ("checkup", "  Health check-up / preventive"), ("other", "  Other")]:
        n = int(cat.get(key, 0))
        add(lab, NA, NA, f"{n} ({100 * n / len(B):.1f})")

    # --- Primary outcome ---
    def pos_cell(d, col):
        return f"{int(d[col].sum())} ({100 * d[col].mean():.1f})"

    tabA = np.array([[pwh["pt_pos"].sum(), len(pwh) - pwh["pt_pos"].sum()],
                     [A["pt_pos"].sum(), len(A) - A["pt_pos"].sum()]])
    chi2_a, p_pos_a = stats.chi2_contingency(tabA, correction=False)[:2]
    tabB = np.array([[pwh["pt_pos"].sum(), len(pwh) - pwh["pt_pos"].sum()],
                     [B["pt_pos"].sum(), len(B) - B["pt_pos"].sum()]])
    chi2_b, p_pos_b = stats.chi2_contingency(tabB, correction=False)[:2]
    add(f"PT-IgG >={C.CUT_POS} {C.ASSAY_UNIT} - n (%)",
        pos_cell(pwh, "pt_pos"), pos_cell(A, "pt_pos"), pos_cell(B, "pt_pos"),
        fmt_p(p_pos_a), fmt_p(p_pos_b))
    C.check("pwoh_a_pos_n", int(A["pt_pos"].sum()), label="PWoH-A seropositive n")
    C.check_pct("pwoh_a_pos_pct", 100 * A["pt_pos"].mean(), label="PWoH-A seropositive %")
    C.check("table1_a_pos_chi2", round(float(chi2_a), 3), label="PWH vs A chi-square")

    def p_fisher(d):
        t = np.array([[pwh["pt_high"].sum(), len(pwh) - pwh["pt_high"].sum()],
                      [d["pt_high"].sum(), len(d) - d["pt_high"].sum()]])
        return float(stats.fisher_exact(t)[1])
    add(f"PT-IgG >={C.CUT_HIGH} {C.ASSAY_UNIT} - n (%)",
        pos_cell(pwh, "pt_high"), pos_cell(A, "pt_high"), pos_cell(B, "pt_high"),
        fmt_p(p_fisher(A)), fmt_p(p_fisher(B)))
    C.check("pwoh_a_high_n", int(A["pt_high"].sum()), label="PWoH-A >=100 n")

    g1, gA, gB = gmc(pwh["pt"]), gmc(A["pt"]), gmc(B["pt"])
    def gcell(g):
        return f"{g[0]:.1f} ({g[1]:.1f}-{g[2]:.1f})"
    add(f"PT-IgG GMC, {C.ASSAY_UNIT} (95% CI)", gcell(g1), gcell(gA), gcell(gB),
        fmt_p(float(stats.ttest_ind(np.log(pwh["pt"]), np.log(A["pt"]), equal_var=False).pvalue)),
        fmt_p(float(stats.ttest_ind(np.log(pwh["pt"]), np.log(B["pt"]), equal_var=False).pvalue)))
    C.check("pwoh_a_gmc", round(gA[0], 1), tol=C.TOL_PCT, label="PWoH-A GMC")

    t1 = pd.DataFrame(rows).fillna("")   # keep empty p-value cells empty in the CSV
    t1.to_csv(C.DIR_NUM / "table1.csv", index=False)
    print(t1.to_string(index=False))
    print(f"\n  [footnote] {C.PWOH_A_UNAVAILABLE_NOTE}")

    # --- Seropositivity by age group: set A is the main figure, set B a sensitivity table ---
    strat = age_strat(A, pwh)
    strat.to_csv(C.DIR_NUM / "table1_agestrat.csv", index=False)
    strat_b = age_strat(B, pwh)
    strat_b.to_csv(C.DIR_NUM / "table1_agestrat_setB.csv", index=False)
    print("\n[By age group, PWoH-A vs PWH]")
    print(strat.to_string(index=False))
    print("\n[By age group, PWoH-B vs PWH (sensitivity)]")
    print(strat_b.to_string(index=False))

    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    x = np.arange(len(strat)); w = 0.36
    ax.bar(x - w / 2, strat["PWoH_pct"], w, label=f"PWoH (n={len(A)})", color=C.PALETTE["pwoh_2025"])
    ax.bar(x + w / 2, strat["PWH_pct"], w, label=f"PWH (n={len(pwh)})", color=C.PALETTE["pwh_2025"])
    for xi, v in zip(x - w / 2, strat["PWoH_pct"]):
        ax.text(xi, v + 1.2, f"{v:.1f}", ha="center", fontsize=8.5)
    for xi, v in zip(x + w / 2, strat["PWH_pct"]):
        ax.text(xi, v + 1.2, f"{v:.1f}", ha="center", fontsize=8.5)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{a}\n({b}/{c})" for a, b, c in
                        zip(strat["age_grp"], strat["PWoH_N"], strat["PWH_N"])], fontsize=8.5)
    ax.set_xlabel("Age group (years) — (PWoH n / PWH n)")
    ax.set_ylabel(f"Anti-PT IgG $\\geq${C.CUT_POS} {C.ASSAY_UNIT} (%)")
    ax.set_title("Seropositivity by age group, PWoH vs PWH (2025)", fontsize=10, loc="left")
    ax.legend(frameon=False, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    stem = C.save_fig(fig, "age_pwoh_vs_pwh")
    plt.close(fig)

    print(f"\n  -> numbers/table1.csv, numbers/table1_agestrat.csv, "
          f"numbers/table1_agestrat_setB.csv, figures/{stem}.pdf")

    return {
        "table1_pwh_pos_pct": round(100 * pwh["pt_pos"].mean(), 1),
        "table1_pwoh_a_pos_pct": round(100 * A["pt_pos"].mean(), 1),
        "table1_pwoh_b_pos_pct": round(100 * B["pt_pos"].mean(), 1),
        "table1_a_pos_chi2": round(float(chi2_a), 3),
        "table1_a_pos_p": float(p_pos_a),
        "table1_b_pos_chi2": round(float(chi2_b), 3),
        "table1_b_pos_p": float(p_pos_b),
        "table1_pwh_gmc": round(g1[0], 1),
        "table1_pwoh_a_gmc": round(gA[0], 1),
        "table1_pwoh_b_gmc": round(gB[0], 1),
        "table1_pwoh_sets": {"A": len(A), "B": len(B)},
    }


if __name__ == "__main__":
    main()
    sys.exit(0)
