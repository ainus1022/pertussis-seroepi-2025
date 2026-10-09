"""
08_biobank_descriptive.py - description of the serum-bank component of PWoH set A.

The 149 serum-bank samples (healthy volunteers) are part of PWoH set A. This
script describes that component on its own, by age decade, so that the two
sources of set A can be reported separately.

Input : output/frozen_biobank.csv   (private)
Output: numbers/biobank_n149.csv
"""
import sys
import numpy as np
import pandas as pd

import config as C

# Age-decade labels in the source workbook -> English
AGE_GROUP_EN = {"20代": "20-29", "30代": "30-39", "40代": "40-49",
                "50代": "50-59", "60代以上": "60+"}


def gm(s):
    return float(np.exp(np.log(s.astype(float)).mean()))


def main():
    print("=" * 72)
    print("08_biobank_descriptive - serum-bank PWoH (healthy volunteers), n=149")
    print("=" * 72)

    bb = pd.read_csv(C.FROZEN_BIOBANK)
    print(f"  dispensed {C.BIOBANK['n_dispensed']} -> measured {len(bb)}")
    print(f"  at the upper limit of quantification (>= {C.ASSAY_ULOQ}): "
          f"{int((bb['pt'] >= C.ASSAY_ULOQ).sum())} (analysed as {C.ASSAY_ULOQ})")
    print(f"  OK  seropositive {int(bb['pt_pos'].sum())}/{len(bb)} ({100 * bb['pt_pos'].mean():.1f}%), "
          f"high titre {int(bb['pt_high'].sum())}, GMC {gm(bb['pt']):.1f}")

    summ = (bb.groupby("age_grp_ja")
            .agg(n=("pt", "size"), male=("male", "sum"), pos=("pt_pos", "sum"),
                 pos_pct=("pt_pos", lambda s: round(100 * s.mean(), 1)),
                 gmc=("pt", lambda s: round(gm(s), 1)))
            .reset_index())
    summ.insert(0, "age_grp", summ.pop("age_grp_ja").map(AGE_GROUP_EN))
    if summ["age_grp"].isna().any():
        raise AssertionError("unexpected age-decade label in the serum-bank file")
    summ.to_csv(C.DIR_NUM / "biobank_n149.csv", index=False)
    print("\n[By age decade]")
    print(summ.to_string(index=False))
    print("\n  -> numbers/biobank_n149.csv")

    return {"biobank_n": len(bb),
            "biobank_pos_pct": round(100 * bb["pt_pos"].mean(), 1),
            "biobank_gmc": round(gm(bb["pt"]), 1)}


if __name__ == "__main__":
    main()
    sys.exit(0)
