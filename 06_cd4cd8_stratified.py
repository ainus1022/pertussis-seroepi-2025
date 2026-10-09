"""
06_cd4cd8_stratified.py - PWH: baseline (2021) CD4/CD8 ratio and 2025 serostatus.

Strata: ratio <= 0.5 versus > 0.5 (a ratio of exactly 0.5 is in the low stratum).
CD8 counts were not measured in 2025, so the ratio is the 2021 baseline value.

Input : output/analysis_master.csv   (private)
Output: numbers/cd4cd8.csv, figures/Fig5_*
"""
import sys
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import config as C


def gmc(s):
    x = np.log(s.astype(float))
    n = len(x)
    se = x.std(ddof=1) / np.sqrt(n)
    t = stats.t.ppf(0.975, n - 1)
    return np.exp(x.mean()), np.exp(x.mean() - t * se), np.exp(x.mean() + t * se)


def main():
    print("=" * 72)
    print("06_cd4cd8_stratified - baseline CD4/CD8 ratio strata")
    print("=" * 72)

    m = pd.read_csv(C.MASTER)
    p = m[(m["hiv"] == 1) & (m["cd48_linked"] == 1)].copy()
    C.check("cd48_pairs", len(p), label="CD4/CD8 analysis n")

    p["stratum"] = np.where(p["cd48_ratio_2021"] <= C.CD48_CUT, "low", "high")
    lo = p[p["stratum"] == "low"]
    hi = p[p["stratum"] == "high"]

    n_half = int((p["cd48_ratio_2021"] == C.CD48_CUT).sum())
    C.check("cd48_exactly_half_n", n_half, label="ratio exactly at the cut-off n")
    half_pos = int(p.loc[p["cd48_ratio_2021"] == C.CD48_CUT, "pt_pos"].sum())
    print(f"  ratio exactly {C.CD48_CUT}: {n_half} (seropositive {half_pos}) -> low stratum")

    C.check("cd48_low_n", len(lo), label="low stratum N")
    C.check("cd48_high_n", len(hi), label="high stratum N")
    C.check("cd48_low_pos", int(lo["pt_pos"].sum()), label="low stratum seropositive n")
    C.check("cd48_high_pos", int(hi["pt_pos"].sum()), label="high stratum seropositive n")
    C.check_pct("cd48_low_pos_pct", 100 * lo["pt_pos"].mean(), label="low stratum seropositive %")
    C.check_pct("cd48_high_pos_pct", 100 * hi["pt_pos"].mean(), label="high stratum seropositive %")
    C.check("cd48_low_median", lo["pt"].median(), label="low stratum median titre")
    C.check("cd48_high_median", hi["pt"].median(), label="high stratum median titre")

    tab = np.array([[lo["pt_pos"].sum(), len(lo) - lo["pt_pos"].sum()],
                    [hi["pt_pos"].sum(), len(hi) - hi["pt_pos"].sum()]])
    or_, p_fisher = stats.fisher_exact(tab)
    chi2, p_chi2 = stats.chi2_contingency(tab, correction=False)[:2]
    chi2_cc, p_chi2_cc = stats.chi2_contingency(tab, correction=True)[:2]   # Yates
    # SciPy's mannwhitneyu applies a continuity correction by default; Stata's
    # ranksum does not. Small differences in p come from this default.
    p_mwu = stats.mannwhitneyu(lo["pt"], hi["pt"], alternative="two-sided").pvalue
    C.check("cd48_fisher_p", float(p_fisher), label="Fisher p")
    C.check("cd48_mwu_p", float(p_mwu), label="Mann-Whitney p")

    g_lo = gmc(lo["pt"]); g_hi = gmc(hi["pt"])
    out = pd.DataFrame([
        {"stratum": f"CD4/CD8 <={C.CD48_CUT}", "N": len(lo),
         "seropositive_n": int(lo["pt_pos"].sum()),
         "seropositive_pct": round(100 * lo["pt_pos"].mean(), 1),
         "GMC": round(g_lo[0], 1), "GMC_lci": round(g_lo[1], 1), "GMC_uci": round(g_lo[2], 1),
         "PT_median": lo["pt"].median()},
        {"stratum": f"CD4/CD8 >{C.CD48_CUT}", "N": len(hi),
         "seropositive_n": int(hi["pt_pos"].sum()),
         "seropositive_pct": round(100 * hi["pt_pos"].mean(), 1),
         "GMC": round(g_hi[0], 1), "GMC_lci": round(g_hi[1], 1), "GMC_uci": round(g_hi[2], 1),
         "PT_median": hi["pt"].median()},
    ])
    out["OR_fisher"] = [round(float(or_), 3), ""]
    out["P_fisher"] = [round(float(p_fisher), 4), ""]
    out["chi2"] = [round(float(chi2), 4), ""]
    out["P_chi2"] = [round(float(p_chi2), 4), ""]
    out["P_chi2_yates"] = [round(float(p_chi2_cc), 4), ""]
    out["P_mwu"] = [round(float(p_mwu), 4), ""]
    out.to_csv(C.DIR_NUM / "cd4cd8.csv", index=False)
    print("\n" + out.to_string(index=False))
    print(f"\n  OK  Fisher p={p_fisher:.3f} / chi2={chi2:.3f} (p={p_chi2:.3f}, "
          f"Yates p={p_chi2_cc:.3f}) / Mann-Whitney p={p_mwu:.3f}")

    # --- Figure ---
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.0))
    ax = axes[0]
    vals = [100 * lo["pt_pos"].mean(), 100 * hi["pt_pos"].mean()]
    labs = [f"$\\leq${C.CD48_CUT}\n(n={len(lo)})", f">{C.CD48_CUT}\n(n={len(hi)})"]
    ax.bar(labs, vals, color=[C.PALETTE["pwh_2025"], C.PALETTE["pwh_2025"]], alpha=None, width=0.55)
    ax.patches[0].set_alpha(0.45)  # same population colour, lighter tint for the low stratum
    for i, v in enumerate(vals):
        ax.text(i, v + 1.2, f"{v:.1f}%", ha="center", fontsize=10)
    ax.set_ylabel(f"Anti-PT IgG $\\geq${C.CUT_POS} {C.ASSAY_UNIT} in 2025 (%)")
    ax.set_xlabel("Baseline CD4/CD8 ratio (2021)")
    ax.set_ylim(0, 70)
    ax.set_title(f"A. Seropositivity (Fisher p={p_fisher:.3f})", fontsize=10, loc="left")
    ax.spines[["top", "right"]].set_visible(False)

    ax = axes[1]
    bp = ax.boxplot([lo["pt"], hi["pt"]], tick_labels=labs, widths=0.5,
                    showfliers=False, patch_artist=True)
    for patch, a in zip(bp["boxes"], [0.30, 0.60]):
        patch.set_facecolor(C.PALETTE["pwh_2025"]); patch.set_alpha(a)
    for med in bp["medians"]:
        med.set_color(C.PALETTE["ink"])
    rng = np.random.default_rng(0)
    for i, s in enumerate([lo["pt"], hi["pt"]], start=1):
        ax.scatter(rng.normal(i, 0.055, len(s)), s, s=7, alpha=0.35, color=C.PALETTE["ink"], lw=0)
    ax.set_yscale("log")
    ax.axhline(C.CUT_POS, color=C.PALETTE["ref_line"], ls="--", lw=1)
    ax.set_ylabel(f"Anti-PT IgG in 2025 ({C.ASSAY_UNIT}, log scale)")
    ax.set_xlabel("Baseline CD4/CD8 ratio (2021)")
    ax.set_title(f"B. Titre distribution (Mann-Whitney p={p_mwu:.3f})", fontsize=10, loc="left")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    stem = C.save_fig(fig, "cd48")
    plt.close(fig)
    print(f"  -> numbers/cd4cd8.csv, figures/{stem}.pdf")

    return {"cd48_pairs": len(p), "cd48_low_pct": round(100 * lo["pt_pos"].mean(), 1),
            "cd48_high_pct": round(100 * hi["pt_pos"].mean(), 1),
            "cd48_fisher_p": round(float(p_fisher), 4),
            "cd48_mwu_p": round(float(p_mwu), 4)}


if __name__ == "__main__":
    main()
    sys.exit(0)
