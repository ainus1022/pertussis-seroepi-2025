"""
14_agestrat_thresholds.py - anti-PT IgG >= 10 / >= 50 / >= 100 EU/mL by age group

Proportion above each threshold by 10-year age group, in PWoH set A and in PWH.
Reading: >= 10 rises with age (lifetime cumulative exposure and whole-cell priming),
whereas >= 50 is nearly flat across ages (recent infection during the epidemic).
The contrast itself is internal evidence that >= 50 marks recent infection.
At >= 100 each stratum has single-digit events, so it is described only.

Age bands are the same as in 03_analysis1_jihs_vs_pwoh.py (20-29 ... 70+), so the
>= 10 gradient can be read next to analysis1_agestrat.csv; the PWoH >= 10 counts
are checked against that file. Thresholds: config.PRIMING_THRESHOLD_CUTS.

Output: numbers/agestrat_thresholds.csv
"""
import sys
import numpy as np
import pandas as pd

import config as C
import lib

# same bands as 03
AGE_BINS = [20, 30, 40, 50, 60, 70, 200]
AGE_LABELS = ["20-29", "30-39", "40-49", "50-59", "60-69", "70+"]


def _binned(df):
    df = df.copy()
    df["age_grp"] = pd.cut(df["age"], bins=AGE_BINS, labels=AGE_LABELS, right=False)
    return df


def _rows(df, population):
    out = []
    for g in AGE_LABELS:
        s = df[df["age_grp"] == g]
        if not len(s):
            continue
        r = {"population": population, "age_grp": g, "N": len(s),
             "age_median": float(s["age"].median())}
        for c in C.PRIMING_THRESHOLD_CUTS:
            r[f"ge{c}_n"] = int((s["pt"] >= c).sum())
            r[f"ge{c}_pct"] = round(100 * float((s["pt"] >= c).mean()), 1)
        out.append(r)
    return out


def _spread(t, population, cut):
    """Range across age bands for one threshold (flatness)."""
    v = t.loc[t["population"] == population, f"ge{cut}_pct"].astype(float)
    return float(v.min()), float(v.max())


def main():
    print("=" * 72)
    print("14_agestrat_thresholds - thresholds by age group")
    print("=" * 72)

    pwoh = _binned(lib.load_pooled_pwoh("14_agestrat_thresholds"))
    assert len(pwoh) == C.EXPECTED["n_pooled_pwoh_2025"], len(pwoh)

    m = pd.read_csv(C.MASTER, parse_dates=["dob"])
    pwh = m[m["hiv"] == 1].copy()
    assert len(pwh) == C.EXPECTED["n_pwh_2025"], len(pwh)
    pwh = _binned(pwh)

    t = pd.DataFrame(_rows(pwoh, "PWoH set A") + _rows(pwh, "PWH"))

    # --- the PWoH >= 10 counts must equal numbers/analysis1_agestrat.csv
    ref = pd.read_csv(C.DIR_NUM / "analysis1_agestrat.csv")
    got = t[t["population"] == "PWoH set A"].set_index("age_grp")
    for _, r in ref.iterrows():
        g = str(r["age_grp"])
        if g not in got.index:
            raise AssertionError(f"age band {g} missing from the new table")
        for a, b in [("PWoH_N", "N"), ("PWoH_pos", "ge10_n")]:
            if int(got.loc[g, b]) != int(r[a]):
                raise AssertionError(
                    f"{g}: {b}={int(got.loc[g, b])} != analysis1_agestrat {a}={int(r[a])}")
    print("  OK: PWoH >=10 by age matches numbers/analysis1_agestrat.csv")

    # --- check >= 50 / >= 100 against config.EXPECTED
    key_age = {"20-29": "20s", "30-39": "30s", "40-49": "40s", "50-59": "50s", "60-69": "60s", "70+": "70p"}
    for pop, tag in [("PWoH set A", "pwoh"), ("PWH", "pwh")]:
        for _, r in t[t["population"] == pop].iterrows():
            a = key_age[str(r["age_grp"])]
            for c in (50, 100):
                C.check(f"agethr_{tag}_{a}_ge{c}_n", int(r[f"ge{c}_n"]), label=f"{pop} {r['age_grp']} >={c} n")
                C.check_pct(f"agethr_{tag}_{a}_ge{c}_pct", float(r[f"ge{c}_pct"]),
                            label=f"{pop} {r['age_grp']} >={c} %")

    t.to_csv(C.DIR_NUM / "agestrat_thresholds.csv", index=False)
    print("\n" + t.to_string(index=False))

    print("\n[range across age bands]")
    for pop in ["PWoH set A", "PWH"]:
        for c in C.PRIMING_THRESHOLD_CUTS:
            lo, hi = _spread(t, pop, c)
            print(f"  {pop:11s} >={c:3d}: {lo:.1f}-{hi:.1f}%  (width {hi - lo:.1f} points)")
    print("\n  ! >=100: single-digit events per band -> described only")
    print(f"\n  -> numbers/agestrat_thresholds.csv")
    return {"agestrat_thresholds": t.to_dict("records")}


if __name__ == "__main__":
    main()
    sys.exit(0)
