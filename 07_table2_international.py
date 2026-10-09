"""
07_table2_international.py - Table 2: comparison with published studies (PWoH / PWH).

Published values are held in config.TABLE2_* and cited by PMID; only the rows
for the present study are computed from the data.

The PWoH row of the present study uses set A (n = 450); set B (n = 301) is
shown as a sensitivity row.

Input : output/analysis_master.csv, output/frozen_biobank.csv, output/frozen_2021.csv (private)
Output: numbers/table2A_pwoh.csv, numbers/table2B_pwh.csv, numbers/table2_footnotes.txt
"""
import sys
import numpy as np
import pandas as pd

import config as C
import lib

CALLER = "07_table2_international"

# Units differ by row: published studies in IU/mL (as published), this study in
# EU/mL. No conversion is applied (config section 2).
COLS_A = ["Study", "Region (adult Tdap)", "n", "Age (y)", "Cutoff", "Unit",
          "Seropositive %", "GMC", "PMID"]
COLS_B = ["Study", "Region", "n", "Age", "CD4", "Cutoff", "Unit",
          "Seropositive %", "GMC", "PMID"]


def with_unit(rows, unit, at):
    """Insert the unit into each literature row at position `at`."""
    return [r[:at] + (unit,) + r[at:] for r in rows]


def gm(s):
    return float(np.exp(np.log(s.astype(float)).mean()))


def main():
    print("=" * 72)
    print("07_table2_international - Table 2")
    print("=" * 72)

    m = pd.read_csv(C.MASTER)
    pwh = m[m["hiv"] == 1]
    pwoh = lib.load_pwoh_set(CALLER)              # A (450)
    pwoh_b = lib.load_pwoh_set(CALLER, "B")       # B (301), sensitivity row
    print(f"  PWoH = set {pwoh.attrs['pwoh_set']} (n={len(pwoh)}); "
          f"set B (n={len(pwoh_b)}) as a sensitivity row")

    # --- Table 2A: PWoH ---
    cuts = "/".join(str(c) for c in C.CUT_MULTI)

    def pwoh_row(d, label):
        return (label, "Japan (no adult Tdap)", str(len(d)),
                f"{int(d['age'].min())}-{int(d['age'].max())}", cuts, C.ASSAY_UNIT,
                " / ".join(f"{100 * (d['pt'] >= c).mean():.1f}" for c in C.CUT_MULTI),
                f"{gm(d['pt']):.1f}", "-")

    row_a = pwoh_row(pwoh, "This study 2025 (PWoH, set A: hospital + serum bank)")
    row_a_b = pwoh_row(pwoh_b, "  (this study PWoH, set B: hospital only)")
    t2a = pd.DataFrame(with_unit(C.TABLE2_PWOH, C.LIT_UNIT, 5) + [row_a, row_a_b],
                       columns=COLS_A)

    # --- Table 2B: PWH (also at the cut-offs used by the comparator studies) ---
    row_b = ("This study 2025 (PWH)", "Japan", str(len(pwh)),
             f"{pwh['age'].median():.0f} (median)", f"{pwh['cd4'].median():.0f} (median)",
             str(C.CUT_POS), C.ASSAY_UNIT, f"{100 * pwh['pt_pos'].mean():.1f}",
             f"{gm(pwh['pt']):.1f}", "-")
    hi_lab = " / ".join(str(c) for c in C.TABLE2B_HIGH_CUTS)
    row_b_hi = ("  (this study PWH, high cutoffs)", "Japan", str(len(pwh)), "", "",
                hi_lab, C.ASSAY_UNIT,
                " / ".join(f"{100 * (pwh['pt'] >= c).mean():.1f}" for c in C.TABLE2B_HIGH_CUTS),
                f"{gm(pwh['pt']):.1f}", "-")
    # The 2021 row is computed from the frozen 2021 file (all 471 with a titre);
    # the analysis file would restrict it to those also sampled in 2025.
    f21 = pd.read_csv(C.FROZEN_2021)
    p21 = f21[f21["in_pt_cohort"] == 1]
    C.check("n_pwh_2021", len(p21), label="Table 2B, 2021 row n")
    row_b_21 = ("  (this study PWH, 2021)", "Japan", str(len(p21)), "", "",
                hi_lab, C.ASSAY_UNIT,
                " / ".join(f"{100 * (p21['pt'] >= c).mean():.1f}" for c in C.TABLE2B_HIGH_CUTS),
                f"{gm(p21['pt']):.1f}", "-")
    t2b = pd.DataFrame(with_unit(C.TABLE2_PWH, C.LIT_UNIT, 5) + [row_b, row_b_hi, row_b_21],
                       columns=COLS_B)

    t2a.to_csv(C.DIR_NUM / "table2A_pwoh.csv", index=False)
    t2b.to_csv(C.DIR_NUM / "table2B_pwh.csv", index=False)
    print("\n[Table 2A] PT-IgG seroprevalence in adults WITHOUT HIV")
    print(t2a.to_string(index=False))
    print("\n[Table 2B] PT-IgG seroprevalence in people WITH HIV")
    print(t2b.to_string(index=False))

    (C.DIR_NUM / "table2_footnotes.txt").write_text(
        "Table 2 footnotes\n"
        "-----------------\n"
        f"[units] {C.TABLE2_UNIT_FOOTNOTE}\n"
        "[Greve-Isdahl 2025] Pregnant women and mother-infant pairs, not a general adult population.\n"
        "[Wehlin 2021] 14 European countries with heterogeneous adult Tdap policies and coverage.\n"
        "[Moosa 2025] Household cohorts (PHIRST); not HIV-specific.\n"
        "[Troy 2016] Recent infection was defined as anti-PT IgG >=94 IU/mL "
        "(cut-off validated by Baughman et al., PMID 15539504). The reported 4.3% is "
        "restricted to participants who had not received Tdap in the preceding six months; "
        "the unrestricted figure was 15/299 (5.0%). Adult Tdap is not used in Japan, so no "
        "such restriction was applied in the present study.\n", encoding="utf-8")
    print("\n  -> numbers/table2A_pwoh.csv, numbers/table2B_pwh.csv, numbers/table2_footnotes.txt")

    return {
        "table2_pwoh_set": pwoh.attrs["pwoh_set"],
        "table2_pwoh_multicut": {str(c): round(100 * (pwoh["pt"] >= c).mean(), 1)
                                 for c in C.CUT_MULTI},
        "table2_pwoh_multicut_setB": {str(c): round(100 * (pwoh_b["pt"] >= c).mean(), 1)
                                      for c in C.CUT_MULTI},
        "table2_pwoh_gmc": round(gm(pwoh["pt"]), 1),
        "table2_pwoh_gmc_setB": round(gm(pwoh_b["pt"]), 1),
        "table2_pwh_gmc": round(gm(pwh["pt"]), 1),
    }


if __name__ == "__main__":
    main()
    sys.exit(0)
