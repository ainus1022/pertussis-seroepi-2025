"""
symptoms_descriptive.py - respiratory symptoms at blood draw and anti-PT IgG
(descriptive; supports the Limitations).

Symptoms come from a retrospective chart review and were not elicited
systematically: no note means "not recorded", not "no symptoms". The cough
category was hand-coded after the antibody results were known. Nothing here
enters the primary analyses.

Reported:
  - participants with any recorded respiratory symptom vs none
  - the subset with cough or sputum recorded
  - how many high-titre participants had a recorded symptom
  - the HIV odds ratio of Analysis 2 with a symptom term added

Input : output/analysis_master.csv   (private)
Output: numbers/symptoms_descriptive.csv, numbers/symptoms_model_check.csv
"""
import sys
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm

import config as C

CUT_MID = 40   # intermediate cut-off used for description


def gmc(x):
    return float(np.exp(np.log(np.clip(x, 1, None)).mean()))


def gmc_ci(x):
    l = np.log(np.clip(x, 1, None))
    se = l.std(ddof=1) / np.sqrt(len(l))
    t = stats.t.ppf(0.975, len(l) - 1)
    return float(np.exp(l.mean() - t * se)), float(np.exp(l.mean() + t * se))


def describe(label, s, ref=None):
    lo, hi = gmc_ci(s["pt"])
    row = {
        "group": label, "n": len(s),
        "GMC": round(gmc(s["pt"]), 1), "GMC_lci": round(lo, 1), "GMC_uci": round(hi, 1),
        "median": float(s["pt"].median()),
        f"ge{C.CUT_POS}_n": int(s["pt_pos"].sum()),
        f"ge{C.CUT_POS}_pct": round(100 * s["pt_pos"].mean(), 1),
        f"ge{CUT_MID}_pct": round(100 * (s["pt"] >= CUT_MID).mean(), 1),
        f"ge{C.CUT_HIGH}_n": int((s["pt"] >= C.CUT_HIGH).sum()),
        "P_mwu_vs_none": "", "P_fisher_seropositive_vs_none": "",
    }
    if ref is not None:
        row["P_mwu_vs_none"] = round(float(stats.mannwhitneyu(s["pt"], ref["pt"]).pvalue), 4)
        tab = [[s["pt_pos"].sum(), len(s) - s["pt_pos"].sum()],
               [ref["pt_pos"].sum(), len(ref) - ref["pt_pos"].sum()]]
        row["P_fisher_seropositive_vs_none"] = round(float(stats.fisher_exact(tab)[1]), 4)
    return row


def main():
    print("=" * 72)
    print("symptoms_descriptive - respiratory symptoms at blood draw (descriptive)")
    print("=" * 72)

    d = pd.read_csv(C.MASTER)
    C.check("n_symptom_linked", len(d), label="2025 records with a symptom record")
    C.check("n_any_sx", int(d["any_sx"].sum()), label="any respiratory symptom recorded n")
    C.check("n_sx_cough", int(d["sx_cough"].sum()), label="cough or sputum recorded n")

    ref = d[d["any_sx"] == 0]
    rows = [describe("No symptom recorded (reference)", ref),
            describe("Any respiratory symptom recorded", d[d["any_sx"] == 1], ref),
            describe("Cough or sputum recorded", d[d["sx_cough"] == 1], ref)]
    out = pd.DataFrame(rows)
    out.to_csv(C.DIR_NUM / "symptoms_descriptive.csv", index=False)
    print(out.to_string(index=False))

    pwh, pwoh = d[d.hiv == 1], d[d.hiv == 0]
    p_grp = stats.fisher_exact(pd.crosstab(d.hiv, d.any_sx).values)[1]
    print(f"\n  any symptom recorded: PWH {int(pwh['any_sx'].sum())}/{len(pwh)} "
          f"({100 * pwh['any_sx'].mean():.1f}%) vs PWoH {int(pwoh['any_sx'].sum())}/{len(pwoh)} "
          f"({100 * pwoh['any_sx'].mean():.1f}%), Fisher p={p_grp:.2f}")
    for cut in (CUT_MID, C.CUT_HIGH):
        h = d[d["pt"] >= cut]
        print(f"  anti-PT IgG >= {cut}: n={len(h)}, with a recorded symptom {int(h['any_sx'].sum())} "
              f"({100 * h['any_sx'].mean():.0f}%)")

    # --- Does adding a symptom term change the HIV odds ratio of Analysis 2? ---
    cc = d[d["cc_analysis2"] == 1]
    mrows = []
    for extra, label in [([], "Analysis 2 model"), (["any_sx"], "+ any symptom"),
                         (["sx_cough"], "+ cough or sputum")]:
        cols = {"hiv": cc["hiv"], "age10": cc["age"] / C.SCALE_AGE, "male": cc["male"],
                "lym500": cc["lym"] / C.SCALE_LYM}
        for e in extra:
            cols[e] = cc[e]
        r = sm.Logit(cc["pt_pos"], sm.add_constant(pd.DataFrame(cols))).fit(disp=0)
        lo_, hi_ = np.exp(r.conf_int().loc["hiv"])
        mrows.append({"model": label, "N": int(r.nobs),
                      "HIV_aOR": round(float(np.exp(r.params["hiv"])), 3),
                      "LCI": round(float(lo_), 3), "UCI": round(float(hi_), 3),
                      "P": round(float(r.pvalues["hiv"]), 4)})
    mc = pd.DataFrame(mrows)
    mc.to_csv(C.DIR_NUM / "symptoms_model_check.csv", index=False)
    print("\n[HIV odds ratio with a symptom term added]")
    print(mc.to_string(index=False))
    print("\n  -> numbers/symptoms_descriptive.csv, numbers/symptoms_model_check.csv")


if __name__ == "__main__":
    main()
    sys.exit(0)
