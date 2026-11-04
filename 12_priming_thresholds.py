"""
12_priming_thresholds.py - Supplementary Table S1: priming-vaccine generation x
anti-PT IgG threshold

Proportion with anti-PT IgG >= 10 / >= 50 / >= 100 EU/mL by priming-vaccine
generation (wP born <= 1974, transition 1975-1981, aP born >= 1982), in PWoH set A
and in PWH 2025. Thresholds: config.PRIMING_THRESHOLD_CUTS.

Reading: if the proportion >= 50 is similar across generations during the epidemic,
exposure was comparable, and the difference at >= 10 reflects the antibody that
persists after exposure rather than exposure itself. At >= 100 the aP generation
has a single event in each population, so that comparison has no power.

aP vs wP is tested with Fisher's exact test (transition cohort excluded).
The >= 10 results are already checked against the keys of 09 and 10; only >= 50
and >= 100 are checked here.

Output: numbers/priming_thresholds.csv
        <private>/output/priming_thresholds_cc.csv  (population, generation, pt;
        input for verify_in_stata.do section 7)
"""
import sys
import numpy as np
import pandas as pd
from scipy import stats

import config as C
import lib

ORDER = ["wP", "transition", "aP"]


def _rows(df, population):
    out = []
    for g in ORDER:
        s = df[df["generation"] == g]
        if not len(s):
            continue
        r = {"population": population, "generation": g, "N": len(s),
             "age_median": float(s["age"].median())}
        for c in C.PRIMING_THRESHOLD_CUTS:
            r[f"ge{c}_n"] = int((s["pt"] >= c).sum())
            r[f"ge{c}_pct"] = round(100 * float((s["pt"] >= c).mean()), 1)
        out.append(r)
    return out


def _tests(df):
    a = df[df["generation"] == "aP"]
    w = df[df["generation"] == "wP"]
    res = {}
    for c in C.PRIMING_THRESHOLD_CUTS:
        ct = np.array([[(a["pt"] >= c).sum(), (a["pt"] < c).sum()],
                       [(w["pt"] >= c).sum(), (w["pt"] < c).sum()]])
        orv, p = stats.fisher_exact(ct)
        res[c] = (float(orv), float(p), ct.tolist())
    return res


def main():
    print("=" * 72)
    print("12_priming_thresholds - Suppl. Table S1: generation x anti-PT IgG threshold")
    print("=" * 72)

    pwoh = lib.load_pooled_pwoh("12_priming_thresholds")
    assert len(pwoh) == C.EXPECTED["n_pwoh_2025"] + C.BIOBANK["n_analyzed"], len(pwoh)

    m = pd.read_csv(C.MASTER, parse_dates=["dob"])
    pwh = m[m["hiv"] == 1].copy()
    assert len(pwh) == C.EXPECTED["n_pwh_2025"], len(pwh)
    pwh["generation"] = pwh["dob"].dt.year.map(lib.vaccine_era)

    rows, tests = [], {}
    for df, pop, tag in [(pwoh, "PWoH set A", "pwoh"), (pwh, "PWH", "pwh")]:
        rows += _rows(df, pop)
        tests[tag] = _tests(df)

    t = pd.DataFrame(rows)
    # Fisher results (aP vs wP, transition excluded) on the aP row of each population
    for tag, pop in [("pwoh", "PWoH set A"), ("pwh", "PWH")]:
        for c, (orv, p, ct) in tests[tag].items():
            mask = (t["population"] == pop) & (t["generation"] == "aP")
            t.loc[mask, f"ge{c}_OR_aP_vs_wP"] = round(orv, 3)
            t.loc[mask, f"ge{c}_P_fisher"] = round(p, 4)

    t.to_csv(C.DIR_NUM / "priming_thresholds.csv", index=False)
    # Stata input: three columns only, no identifiers
    cc = pd.concat([
        pwoh.assign(population="PWoH")[["population", "generation", "pt"]],
        pwh.assign(population="PWH")[["population", "generation", "pt"]],
    ], ignore_index=True)
    cc.to_csv(C.DIR_OUT / "priming_thresholds_cc.csv", index=False)

    print("\n" + t.to_string(index=False))
    for tag, pop in [("pwoh", "PWoH set A"), ("pwh", "PWH")]:
        print(f"\n  {pop}: aP vs wP")
        for c, (orv, p, ct) in tests[tag].items():
            print(f"    >={c:3d}: OR={orv:.3f}  Fisher p={p:.4f}  {ct}")

    # ---- check against config.EXPECTED
    def pick(pop, gen, col):
        return float(t.loc[(t["population"] == pop) & (t["generation"] == gen), col].iloc[0])

    for tag, pop in [("pwoh", "PWoH set A"), ("pwh", "PWH")]:
        for gen, gk in [("wP", "wp"), ("aP", "ap"), ("transition", "trans")]:
            for c in C.PRIMING_THRESHOLD_CUTS:
                if c == 10:
                    continue           # >= 10 is checked by 09 / 10
                C.check(f"primthr_{tag}_{gk}_ge{c}_n", int(pick(pop, gen, f"ge{c}_n")),
                        label=f"{pop} {gen} >={c} n")
                C.check_pct(f"primthr_{tag}_{gk}_ge{c}_pct", pick(pop, gen, f"ge{c}_pct"),
                            label=f"{pop} {gen} >={c} %")
        for c in C.PRIMING_THRESHOLD_CUTS:
            if c == 10:
                continue
            C.check(f"primthr_{tag}_ge{c}_fisher_p", tests[tag][c][1],
                    label=f"{pop} aP vs wP >={c} Fisher p")

    print("\n  -> numbers/priming_thresholds.csv")
    return {"priming_thresholds": t.to_dict("records")}


if __name__ == "__main__":
    main()
    sys.exit(0)
