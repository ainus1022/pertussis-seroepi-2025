"""
05_analysis3_mcnemar.py - Analysis 3: PWH, paired sera 2021 -> 2025 (McNemar).

Input : output/analysis_master.csv   (private; pairs linked by study ID + date of birth)
Output: numbers/mcnemar_paired.csv, numbers/mcnemar_agegroup.csv, figures/Fig4_*
"""
import sys
import numpy as np
import pandas as pd
from statsmodels.stats.contingency_tables import mcnemar
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import config as C

AGE_BINS = [0, 40, 50, 60, 200]
AGE_LABELS = ["<40", "40-49", "50-59", "60+"]


def main():
    print("=" * 72)
    print("05_analysis3 - PWH paired 2021 -> 2025, McNemar")
    print("=" * 72)

    m = pd.read_csv(C.MASTER)
    p = m[(m["hiv"] == 1) & (m["paired"] == 1)].copy()
    C.check("mcnemar_pairs", len(p), label="longitudinal pairs n")

    a = p["pt_pos_2021"].astype(int)   # seropositive in 2021
    b = p["pt_pos"].astype(int)        # seropositive in 2025
    n11 = int(((a == 1) & (b == 1)).sum())
    n10 = int(((a == 1) & (b == 0)).sum())   # positive -> negative
    n01 = int(((a == 0) & (b == 1)).sum())   # negative -> positive
    n00 = int(((a == 0) & (b == 0)).sum())

    C.check("mcnemar_pos_2021", int(a.sum()), label="2021 seropositive n")
    C.check_pct("mcnemar_pos_2021_pct", 100 * a.mean(), label="2021 seropositive %")
    C.check("mcnemar_pos_2025", int(b.sum()), label="2025 seropositive n")
    C.check_pct("mcnemar_pos_2025_pct", 100 * b.mean(), label="2025 seropositive %")
    C.check("mcnemar_disc_neg_to_pos", n01, label="negative -> positive")
    C.check("mcnemar_disc_pos_to_neg", n10, label="positive -> negative")

    # The reported statistic uses the continuity correction. Stata's `mcc` reports
    # the uncorrected statistic; this is a difference in defaults, not in method.
    # The exact test agrees between the two packages and is checked alongside.
    r = mcnemar(np.array([[n11, n10], [n01, n00]]), exact=False, correction=True)
    r_nc = mcnemar(np.array([[n11, n10], [n01, n00]]), exact=False, correction=False)
    r_exact = mcnemar(np.array([[n11, n10], [n01, n00]]), exact=True)
    C.check("mcnemar_chi2", float(r.statistic), label="McNemar chi-square (corrected, reported)")
    C.check("mcnemar_p", float(r.pvalue), label="McNemar p (corrected, reported)")
    C.check("mcnemar_p_exact", float(r_exact.pvalue), label="McNemar exact p")
    C.check("mcnemar_chi2_stata_uncorrected", float(r_nc.statistic),
            label="McNemar chi-square (uncorrected, reference)")

    out = pd.DataFrame([{
        "n_pairs": len(p),
        "pos_2021_n": int(a.sum()), "pos_2021_pct": round(100 * a.mean(), 1),
        "pos_2025_n": int(b.sum()), "pos_2025_pct": round(100 * b.mean(), 1),
        "both_pos": n11, "pos_to_neg": n10, "neg_to_pos": n01, "both_neg": n00,
        "chi2_cc": round(float(r.statistic), 3),          # reported (continuity-corrected)
        "p_chi2_cc": round(float(r.pvalue), 4),           # reported
        "p_exact_binomial": round(float(r_exact.pvalue), 4),
        "chi2_uncorrected": round(float(r_nc.statistic), 3),  # reference (Stata mcc default)
        "p_uncorrected": round(float(r_nc.pvalue), 4),
    }])
    out.to_csv(C.DIR_NUM / "mcnemar_paired.csv", index=False)
    print(out.T.to_string(header=False))
    print(f"\n  OK  {100 * a.mean():.1f}% -> {100 * b.mean():.1f}%  "
          f"chi2={r.statistic:.3f}, p={r.pvalue:.3f} (continuity-corrected)")
    print(f"      reference: uncorrected chi2={r_nc.statistic:.3f}, p={r_nc.pvalue:.4f}; "
          f"exact p={r_exact.pvalue:.4f}")

    # --- By age group (age in 2025) ---
    p["age_grp"] = pd.cut(p["age"], bins=AGE_BINS, labels=AGE_LABELS, right=False)
    rows = []
    for g, s in p.groupby("age_grp", observed=True):
        aa = s["pt_pos_2021"].astype(int); bb = s["pt_pos"].astype(int)
        t = np.array([[int(((aa == 1) & (bb == 1)).sum()), int(((aa == 1) & (bb == 0)).sum())],
                      [int(((aa == 0) & (bb == 1)).sum()), int(((aa == 0) & (bb == 0)).sum())]])
        rr = mcnemar(t, exact=True)
        rows.append({"age_grp": g, "n": len(s),
                     "pos_2021_pct": round(100 * aa.mean(), 1),
                     "pos_2025_pct": round(100 * bb.mean(), 1),
                     "neg_to_pos": int(t[1, 0]), "pos_to_neg": int(t[0, 1]),
                     "p_exact": round(float(rr.pvalue), 4)})
    ag = pd.DataFrame(rows)
    ag.to_csv(C.DIR_NUM / "mcnemar_agegroup.csv", index=False)
    print("\n[By age group]")
    print(ag.to_string(index=False))

    # --- Figure ---
    fig, ax = plt.subplots(figsize=(6.8, 4.0))
    x = np.arange(len(ag)); w = 0.36
    ax.bar(x - w / 2, ag["pos_2021_pct"], w, label="2021", color=C.PALETTE["pwh_2021"])
    ax.bar(x + w / 2, ag["pos_2025_pct"], w, label="2025", color=C.PALETTE["pwh_2025"])
    for xi, v in zip(x - w / 2, ag["pos_2021_pct"]):
        ax.text(xi, v + 1.2, f"{v:.1f}", ha="center", fontsize=8.5)
    for xi, v in zip(x + w / 2, ag["pos_2025_pct"]):
        ax.text(xi, v + 1.2, f"{v:.1f}", ha="center", fontsize=8.5)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{g}\n(n={n})" for g, n in zip(ag["age_grp"], ag["n"])], fontsize=9)
    ax.set_xlabel("Age group in 2025 (years)")
    ax.set_ylabel(f"Anti-PT IgG $\\geq${C.CUT_POS} {C.ASSAY_UNIT} (%)")
    ax.set_title(f"Paired anti-PT IgG seropositivity in PWH, 2021 vs 2025 "
                 f"(n={len(p)} pairs; McNemar p={r.pvalue:.3f})", fontsize=9.5, loc="left")
    ax.legend(frameon=False, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    stem = C.save_fig(fig, "mcnemar_age")
    plt.close(fig)
    print(f"\n  -> numbers/mcnemar_paired.csv, numbers/mcnemar_agegroup.csv, figures/{stem}.pdf")

    return {"mcnemar_pairs": len(p), "mcnemar_chi2": round(float(r.statistic), 3),
            "mcnemar_p": round(float(r.pvalue), 4),
            "mcnemar_pos_2021_pct": round(100 * a.mean(), 1),
            "mcnemar_pos_2025_pct": round(100 * b.mean(), 1)}


if __name__ == "__main__":
    main()
    sys.exit(0)
