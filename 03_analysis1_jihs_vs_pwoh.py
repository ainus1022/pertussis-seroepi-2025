"""
03_analysis1_jihs_vs_pwoh.py - Analysis 1: general population (JIHS 2023) vs PWoH (2025).

Input : data/jihs2023_pt_from_pdf.csv (public), output/analysis_master.csv,
        output/frozen_biobank.csv (private)
Output: numbers/analysis1.csv, numbers/analysis1_agestrat.csv, figures/Fig3_*

The JIHS survey provides age-specific counts by titre bin, not individual
values, so the usable cut-offs are 10, 50 and 100.

The PWoH side is set A (outpatients + serum bank, n = 450): this comparison
needs only anti-PT IgG and age. Because the PWoH sample is older than the
survey sample, the comparison is also shown by age group.
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

CALLER = "03_analysis1_jihs_vs_pwoh"

AGE_BINS = [20, 30, 40, 50, 60, 70, 200]
AGE_LABELS = ["20-29", "30-39", "40-49", "50-59", "60-69", "70+"]


def load_jihs() -> pd.DataFrame:
    """JIHS 2023 survey, Table 3-1, restricted to age >= 20."""
    d = pd.read_csv(C.JIHS_FROM_PDF)
    d = d[d["age_num"] >= 20].copy()
    d["ge10"] = d[["10-49", "50-99", ">=100"]].sum(axis=1)
    d["ge50"] = d[["50-99", ">=100"]].sum(axis=1)
    d["ge100"] = d[">=100"]
    d["age_grp"] = pd.cut(d["age_num"], bins=AGE_BINS, labels=AGE_LABELS, right=False)
    return d


def two_prop(k1, n1, k0, n0):
    """2x2 comparison: chi-square (no continuity correction) and Fisher's exact test."""
    tab = np.array([[k1, n1 - k1], [k0, n0 - k0]])
    chi2, p_chi = stats.chi2_contingency(tab, correction=False)[:2]
    p_fis = stats.fisher_exact(tab)[1]
    return float(chi2), float(p_chi), float(p_fis)


def main():
    print("=" * 72)
    print("03_analysis1 - general population JIHS 2023 vs PWoH 2025")
    print("=" * 72)

    j = load_jihs()

    pwoh = lib.load_pwoh_set(CALLER)
    C.check("analysis1_a_pwoh_n", len(pwoh), label="Analysis 1 PWoH-A n")
    C.check("analysis1_a_pwoh_pos_n", int(pwoh["pt_pos"].sum()), label="Analysis 1 PWoH-A seropositive n")
    print(f"  PWoH = set {pwoh.attrs['pwoh_set']} (n={len(pwoh)}: "
          f"outpatients {int((pwoh['cohort'] == 'main').sum())} + serum bank "
          f"{int((pwoh['cohort'] == 'biobank').sum())})")

    n_j = int(j["Total"].sum())
    C.check("jihs_pdf_n", n_j, label="JIHS n")
    C.check("jihs_pdf_pos_n", int(j["ge10"].sum()), label="JIHS >=10 n")
    C.check_pct("jihs_pdf_pos_pct", 100 * j["ge10"].sum() / n_j, label="JIHS >=10 %")
    C.check("jihs_pdf_high_n", int(j["ge100"].sum()), label="JIHS >=100 n")

    rows = []
    chi2_by_cut = {}
    for cut, jcol in [(10, "ge10"), (50, "ge50"), (100, "ge100")]:
        kj = int(j[jcol].sum())
        kp = int((pwoh["pt"] >= cut).sum())
        chi2, p_chi, p_fis = two_prop(kp, len(pwoh), kj, n_j)
        chi2_by_cut[cut] = round(chi2, 3)
        rows.append({
            "Cutoff (EU/mL)": f">={cut}",
            "JIHS2023 n/N (%)": f"{kj}/{n_j} ({100 * kj / n_j:.1f})",
            f"PWoH2025-{pwoh.attrs['pwoh_set']} n/N (%)":
                f"{kp}/{len(pwoh)} ({100 * kp / len(pwoh):.1f})",
            "Ratio": round((kp / len(pwoh)) / (kj / n_j), 2),
            "chi2": round(chi2, 3),
            "P (chi2)": "<0.001" if p_chi < 0.001 else f"{p_chi:.3f}",
            "P (Fisher)": "<0.001" if p_fis < 0.001 else f"{p_fis:.3f}",
        })
    a1 = pd.DataFrame(rows)
    C.check("analysis1_a_chi2_ge10", chi2_by_cut[10], label="Analysis 1 chi-square >=10")
    C.check("analysis1_a_chi2_ge50", chi2_by_cut[50], label="Analysis 1 chi-square >=50")
    C.check("analysis1_a_chi2_ge100", chi2_by_cut[100], label="Analysis 1 chi-square >=100")

    # --- By age group ---
    pwoh = pwoh.copy()
    pwoh["age_grp"] = pd.cut(pwoh["age"], bins=AGE_BINS, labels=AGE_LABELS, right=False)
    jg = j.groupby("age_grp", observed=True)[["Total", "ge10"]].sum()
    pg = pwoh.groupby("age_grp", observed=True)["pt_pos"].agg(["sum", "count"])
    strat = pd.DataFrame({
        "JIHS_N": jg["Total"], "JIHS_pos": jg["ge10"],
        "JIHS_pct": (100 * jg["ge10"] / jg["Total"]).round(1),
        "PWoH_N": pg["count"], "PWoH_pos": pg["sum"],
        "PWoH_pct": (100 * pg["sum"] / pg["count"]).round(1),
    }).reset_index()
    strat.to_csv(C.DIR_NUM / "analysis1_agestrat.csv", index=False)

    a1.to_csv(C.DIR_NUM / "analysis1.csv", index=False)
    print(a1.to_string(index=False))
    print("\n[By age group]")
    print(strat.to_string(index=False))

    # --- Figure: A, by cut-off; B, by age group ---
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    ax = axes[0]
    cuts = [10, 50, 100]
    jp = [100 * j[c].sum() / n_j for c in ["ge10", "ge50", "ge100"]]
    pp = [100 * (pwoh["pt"] >= c).mean() for c in cuts]
    x = np.arange(len(cuts)); w = 0.36
    ax.bar(x - w / 2, jp, w, label=f"General population, JIHS 2023 (n={n_j})", color=C.PALETTE["jihs_2023"])
    ax.bar(x + w / 2, pp, w, label=f"This study, PWoH 2025 (n={len(pwoh)})", color=C.PALETTE["pwoh_2025"])
    for xi, v in zip(x - w / 2, jp):
        ax.text(xi, v + 1, f"{v:.1f}", ha="center", fontsize=9)
    for xi, v in zip(x + w / 2, pp):
        ax.text(xi, v + 1, f"{v:.1f}", ha="center", fontsize=9)
    ax.set_xticks(x); ax.set_xticklabels([f"$\\geq${c}" for c in cuts])
    ax.set_xlabel(f"Anti-PT IgG cutoff ({C.ASSAY_UNIT})")
    ax.set_ylabel("Seropositive (%)")
    ax.set_title("A. Seropositivity by cutoff")
    ax.legend(fontsize=8, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)

    ax = axes[1]
    x = np.arange(len(strat))
    ax.bar(x - w / 2, strat["JIHS_pct"], w, color=C.PALETTE["jihs_2023"], label="JIHS 2023")
    ax.bar(x + w / 2, strat["PWoH_pct"], w, color=C.PALETTE["pwoh_2025"],
           label=f"PWoH 2025 (n={len(pwoh)})")
    ax.set_xticks(x); ax.set_xticklabels(strat["age_grp"], fontsize=9)
    ax.set_xlabel("Age group (years)")
    ax.set_ylabel(f"Anti-PT IgG $\\geq${C.CUT_POS} {C.ASSAY_UNIT} (%)")
    ax.set_title("B. By age group")
    ax.legend(fontsize=8, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    stem = C.save_fig(fig, "jihs_vs_pwoh")
    plt.close(fig)
    print(f"\n  -> numbers/analysis1.csv, numbers/analysis1_agestrat.csv, figures/{stem}.pdf")

    return {
        "analysis1_jihs_pos_pct": round(100 * j["ge10"].sum() / n_j, 1),
        "analysis1_pwoh_pos_pct": round(100 * pwoh["pt_pos"].mean(), 1),
        "analysis1_jihs_high_pct": round(100 * j["ge100"].sum() / n_j, 2),
        "analysis1_pwoh_high_pct": round(100 * pwoh["pt_high"].mean(), 1),
    }


if __name__ == "__main__":
    main()
    sys.exit(0)
