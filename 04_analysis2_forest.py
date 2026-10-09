"""
04_analysis2_forest.py - Analysis 2: multivariable logistic regression and forest plot.

Main model : outcome = anti-PT IgG >= 10; exposure = HIV; adjusted for age, sex
             and lymphocyte count. Complete cases: 829 - 96 without a lymphocyte
             count = 733.
Sensitivity: 1:1 age matching by the rule fixed in config.MATCH_ALGORITHM.

This is the only PWoH analysis restricted to set B (hospital outpatients). The
model adjusts for lymphocyte count, and serum-bank samples come with sex and age
only, so they cannot enter a complete-case model.

Input : output/analysis_master.csv   (private)
Output: numbers/analysis2.csv, numbers/analysis2_matched.csv, figures/Fig2_*
"""
import sys
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import config as C

CALLER = "04_analysis2_forest"

TERMS = [
    ("hiv", "HIV-positive (PWH vs PWoH)", 1),
    ("age_s", f"Age (per {C.SCALE_AGE} years)", C.SCALE_AGE),
    ("male", "Male sex", 1),
    ("lym_s", f"Lymphocyte count (per {C.SCALE_LYM} /uL)", C.SCALE_LYM),
]


def matched_sensitivity(m: pd.DataFrame, match_sex: bool = False):
    """1:1 matching by config.MATCH_ALGORITHM; conditional OR from discordant pairs."""
    cc = m[m["cc_analysis2"] == 1].copy()
    treat = cc[cc["hiv"] == 1].sort_values(["age", "id"]).to_dict("records")
    ctrl = cc[cc["hiv"] == 0].sort_values(["age", "id"]).to_dict("records")
    used, pairs = set(), []
    for t in treat:
        best, bestd = None, None
        for i, c in enumerate(ctrl):
            if i in used:
                continue
            if match_sex and c["male"] != t["male"]:
                continue
            d = abs(c["age"] - t["age"])
            if d <= C.MATCH_AGE_CALIPER and (bestd is None or d < bestd):
                best, bestd = i, d
        if best is not None:
            used.add(best)
            pairs.append((t["pt_pos"], ctrl[best]["pt_pos"]))
    b = sum(1 for a, c in pairs if a == 1 and c == 0)   # PWH positive, PWoH negative
    d = sum(1 for a, c in pairs if a == 0 and c == 1)   # PWH negative, PWoH positive
    or_m = b / d if d else np.nan
    p = stats.binomtest(b, b + d, 0.5).pvalue if (b + d) else np.nan
    lo, hi = (np.nan, np.nan)
    if b and d:
        se = np.sqrt(1 / b + 1 / d)
        lo, hi = np.exp(np.log(or_m) - 1.96 * se), np.exp(np.log(or_m) + 1.96 * se)
    return {"n_pairs": len(pairs), "discordant_b": b, "discordant_d": d,
            "matched_OR": or_m, "lci": lo, "uci": hi, "p": p}


def main():
    print("=" * 72)
    print("04_analysis2 - multivariable logistic regression and forest plot")
    print("=" * 72)

    which = C.pwoh_set_for(CALLER)
    if which != "B":
        raise AssertionError(
            f"PWoH set for script 04 is '{which}'. The model adjusts for lymphocyte count "
            f"and must be restricted to set B. Check config.PWOH_SET_BY_ANALYSIS.")

    m = pd.read_csv(C.MASTER)
    cc = m[m["cc_analysis2"] == 1].copy()
    C.check("analysis2_n", len(cc), label="Analysis 2 complete cases N")
    dropped = len(m) - len(cc)
    print(f"  PWoH = set B (outpatients, n={C.EXPECTED['n_pwoh_2025']})")
    print(f"  complete cases N={len(cc)} ({dropped} of {len(m)} excluded for a missing covariate)")
    print("    missing: " + ", ".join(
        f"{c} {int(m[c].isna().sum())}" for c in C.COVARIATES))
    print(f"    PWH {int((cc['hiv'] == 1).sum())} / PWoH {int((cc['hiv'] == 0).sum())}")

    cc["age_s"] = cc["age"] / C.SCALE_AGE
    cc["lym_s"] = cc["lym"] / C.SCALE_LYM
    X = sm.add_constant(cc[[t[0] for t in TERMS]].astype(float))
    res = sm.Logit(cc["pt_pos"], X).fit(disp=0)

    ci = res.conf_int()
    rows = []
    for var, label, _ in TERMS:
        rows.append({
            "Covariate": label,
            "aOR": round(float(np.exp(res.params[var])), 3),
            "LCI": round(float(np.exp(ci.loc[var, 0])), 3),
            "UCI": round(float(np.exp(ci.loc[var, 1])), 3),
            "P": round(float(res.pvalues[var]), 4),
        })
    a2 = pd.DataFrame(rows)
    a2.insert(1, "N", [len(cc)] + [""] * (len(TERMS) - 1))

    C.check("analysis2_hiv_aor", a2.loc[0, "aOR"], label="HIV aOR")
    C.check("analysis2_hiv_lci", a2.loc[0, "LCI"], label="HIV aOR lower limit")
    C.check("analysis2_hiv_uci", a2.loc[0, "UCI"], label="HIV aOR upper limit")
    C.check("analysis2_hiv_p", a2.loc[0, "P"], label="HIV aOR p")

    a2.to_csv(C.DIR_NUM / "analysis2.csv", index=False)
    print("\n" + a2.to_string(index=False))

    # --- Sensitivity: age-matched pairs ---
    ms = matched_sensitivity(m, match_sex=False)
    ms_sex = matched_sensitivity(m, match_sex=True)
    C.check("matched_pairs", ms["n_pairs"], label="matched pairs n")
    C.check("matched_or", ms["matched_OR"], label="matched OR")
    C.check("matched_p", ms["p"], label="matched p")

    md = pd.DataFrame([{"variant": "primary (age only)", **ms},
                       {"variant": "variant (age + sex)", **ms_sex}])
    md["algorithm"] = C.MATCH_ALGORITHM
    md.to_csv(C.DIR_NUM / "analysis2_matched.csv", index=False)
    print(f"\n[Sensitivity: 1:1 matching on age within ±{C.MATCH_AGE_CALIPER} years]")
    for lab, r in [("age only (primary)", ms), ("age + sex (variant)", ms_sex)]:
        print(f"  {lab:20s} pairs {r['n_pairs']:3d} / discordant {r['discordant_b']}:{r['discordant_d']} -> "
              f"OR {r['matched_OR']:.2f} ({r['lci']:.2f}-{r['uci']:.2f}) p={r['p']:.3f}")

    # --- Forest plot ---
    fig, ax = plt.subplots(figsize=(7.6, 3.2))
    y = np.arange(len(a2))[::-1]
    ax.errorbar(a2["aOR"], y,
                xerr=[a2["aOR"] - a2["LCI"], a2["UCI"] - a2["aOR"]],
                fmt="o", color=C.PALETTE["ink"], ecolor=C.PALETTE["ink"], capsize=3, markersize=6, lw=1.4)
    ax.axvline(1.0, color=C.PALETTE["ref_line"], ls="--", lw=1)
    ax.set_yticks(y)
    ax.set_yticklabels(a2["Covariate"], fontsize=9)
    ax.set_xscale("log")
    ax.set_xticks([0.25, 0.5, 1, 2])
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.set_xlabel(f"Adjusted odds ratio for anti-PT IgG $\\geq${C.CUT_POS} {C.ASSAY_UNIT} (95% CI)")
    ax.set_xlim(0.2, 2.6)
    for yi, r in zip(y, a2.itertuples()):
        ax.text(2.75, yi, f"{r.aOR:.2f} ({r.LCI:.2f}-{r.UCI:.2f})  p={r.P:.3f}",
                va="center", fontsize=8.5)
    ax.set_title(f"Multivariable logistic regression (complete case, N={len(cc)})",
                 fontsize=10, loc="left")
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    fig.tight_layout()
    stem = C.save_fig(fig, "forest", bbox_inches="tight")
    plt.close(fig)
    print(f"\n  -> numbers/analysis2.csv, numbers/analysis2_matched.csv, figures/{stem}.pdf")

    return {
        "analysis2_pwoh_set": which,
        "analysis2_n": len(cc),
        "analysis2_n_pwh": int((cc["hiv"] == 1).sum()),
        "analysis2_n_pwoh": int((cc["hiv"] == 0).sum()),
        "analysis2_hiv_aor": a2.loc[0, "aOR"],
        "analysis2_hiv_ci": [a2.loc[0, "LCI"], a2.loc[0, "UCI"]],
        "analysis2_hiv_p": a2.loc[0, "P"],
        "analysis2_matched_or": round(float(ms["matched_OR"]), 3),
        "analysis2_matched_p": round(float(ms["p"]), 4),
    }


if __name__ == "__main__":
    main()
    sys.exit(0)
