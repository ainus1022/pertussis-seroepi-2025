"""
01_derive_and_validate.py - derive every flag from the thresholds and build the
single analysis file.

Input : output/frozen_2025.csv, output/frozen_2021.csv   (private)
Output: output/analysis_master.csv                        (private; read by Python and Stata)

All flags are recomputed from CUT_POS / CUT_HIGH. Derived columns in the raw
records are never used.
"""
import sys
import numpy as np
import pandas as pd

import config as C


def classify_dx(dx: str) -> str:
    """Classify the reason for visit of a PWoH outpatient from the diagnosis field.

    Priority: liver > metabolic > preventive > other. Liver disease is the main
    component of this outpatient population and takes precedence when several
    diagnoses are listed. Preventive = health check-up, vaccination, travel
    clinic, infection control.
    """
    s = str(dx)
    if any(k in s for k in C.DX_LIVER):
        return "liver"
    if any(k in s for k in C.DX_METABOLIC):
        return "metabolic"
    if any(k in s for k in C.DX_CHECKUP):
        return "checkup"
    return "other"


def main():
    print("=" * 72)
    print("01_derive_and_validate - flags from thresholds -> analysis file")
    print("=" * 72)

    f25 = pd.read_csv(C.FROZEN_2025, parse_dates=["dob"])
    f21 = pd.read_csv(C.FROZEN_2021, parse_dates=["dob"])

    # --- Flags, from thresholds only ---
    f25["pt_pos"] = (f25["pt"] >= C.CUT_POS).astype(int)
    f25["pt_high"] = (f25["pt"] >= C.CUT_HIGH).astype(int)
    f25["censored"] = (f25["pt"] >= C.ASSAY_ULOQ).astype(int)
    f25["hiv"] = (f25["group"] == "PWH").astype(int)
    for cut in C.CUT_MULTI:
        f25[f"pt_ge{cut}"] = (f25["pt"] >= cut).astype(int)

    # Records without a 2021 titre keep missing flags (not coerced to 0).
    f21["pt_pos"] = np.where(f21["pt"].notna(), (f21["pt"] >= C.CUT_POS).astype(float), np.nan)
    f21["pt_high"] = np.where(f21["pt"].notna(), (f21["pt"] >= C.CUT_HIGH).astype(float), np.nan)

    # --- Respiratory-symptom categories (hand-coded sets in the private configuration).
    #     Descriptive only; not used in the primary analyses.
    sn = f25["sample_no_sx"]
    f25["sx_cough"] = sn.isin(C.private("SX_COUGH")).astype(int)
    f25["sx_prolonged"] = sn.isin(C.private("SX_PROLONGED")).astype(int)
    f25["sx_prolonged_strict"] = sn.isin(C.private("SX_PROLONGED_STRICT")).astype(int)
    f25["sx_altdx"] = sn.isin(C.private("SX_ALTDX")).astype(int)
    for col, key in [("sx_cough", "n_sx_cough"), ("sx_prolonged", "n_sx_prolonged"),
                     ("sx_altdx", "n_sx_altdx")]:
        C.check(key, int(f25[col].sum()), label=f"{col} n")
    if int((f25["sx_cough"] & ~f25["any_sx"].astype(bool)).sum()):
        raise AssertionError("a record coded as cough has no symptom note (sample number error)")

    # --- Reason for visit (PWoH) ---
    f25["dx_cat"] = np.where(f25["hiv"] == 1, "hiv", f25["dx"].map(classify_dx))

    pwh = f25[f25["hiv"] == 1]
    pwoh = f25[f25["hiv"] == 0]

    # --- Analysis sets (re-checked so that this script can run on its own) ---
    C.check("n_pwh_2025", len(pwh), label="PWH 2025")
    C.check("n_pwoh_2025", len(pwoh), label="PWoH 2025")
    C.check("n_pwh_2021", int(f21["in_pt_cohort"].sum()), label="PWH 2021")

    # --- Primary outcome ---
    C.check("pwh_pos_n", int(pwh["pt_pos"].sum()), label="PWH seropositive n")
    C.check_pct("pwh_pos_pct", 100 * pwh["pt_pos"].mean(), label="PWH seropositive %")
    C.check("pwoh_pos_n", int(pwoh["pt_pos"].sum()), label="PWoH seropositive n")
    C.check_pct("pwoh_pos_pct", 100 * pwoh["pt_pos"].mean(), label="PWoH seropositive %")
    C.check("pwh_high_n", int(pwh["pt_high"].sum()), label="PWH >=100 n")
    C.check("pwoh_high_n", int(pwoh["pt_high"].sum()), label="PWoH >=100 n")

    # --- Reason for visit ---
    cat = pwoh["dx_cat"].value_counts()
    C.check("pwoh_liver_n", int(cat.get("liver", 0)), label="PWoH liver disease n")
    C.check("pwoh_metabolic_n", int(cat.get("metabolic", 0)), label="PWoH metabolic disease n")
    C.check("pwoh_checkup_n", int(cat.get("checkup", 0)), label="PWoH preventive visit n")

    print(f"\n  right-censored (anti-PT IgG >= {C.ASSAY_ULOQ}): {int(f25['censored'].sum())} "
          f"(PWH {int(pwh['censored'].sum())} / PWoH {int(pwoh['censored'].sum())})")
    print("  PWoH reason for visit:", {k: int(v) for k, v in cat.items()})
    print(f"  OK  seropositivity PWH {100 * pwh['pt_pos'].mean():.1f}% ({pwh['pt_pos'].sum()}/{len(pwh)}) "
          f"vs PWoH {100 * pwoh['pt_pos'].mean():.1f}% ({pwoh['pt_pos'].sum()}/{len(pwoh)})")

    # --- Longitudinal link (study ID + date of birth) ---
    link = f21[["id", "dob", "pt", "pt_pos", "cd4", "cd8", "cd48_ratio"]].rename(
        columns={"pt": "pt_2021", "pt_pos": "pt_pos_2021", "cd4": "cd4_2021",
                 "cd8": "cd8_2021", "cd48_ratio": "cd48_ratio_2021"})
    master = f25.merge(link, on=["id", "dob"], how="left")
    if len(master) != len(f25):
        raise AssertionError("longitudinal link added rows (ID + date of birth is not 1:1)")

    # paired      : McNemar set (anti-PT IgG in both 2021 and 2025)
    # cd48_linked : CD4/CD8 set (2021 ratio and 2025 titre; a 2021 titre is not
    #               required, so this set is larger than `paired`)
    master["paired"] = master["pt_2021"].notna().astype(int)
    master["cd48_linked"] = (master["cd48_ratio_2021"].notna() & master["pt"].notna()).astype(int)
    C.check("mcnemar_pairs", int(master["paired"].sum()), label="longitudinal pairs n")
    C.check("cd48_pairs", int(master["cd48_linked"].sum()), label="CD4/CD8 analysis n")
    print(f"  OK  longitudinal pairs (2021 x 2025 titre): {int(master['paired'].sum())}")
    print(f"  OK  CD4/CD8 analysis (2021 ratio x 2025 titre): {int(master['cd48_linked'].sum())}")

    # --- Complete-case flag for Analysis 2 (defined in one place) ---
    master["cc_analysis2"] = master[C.COVARIATES].notna().all(axis=1).astype(int)
    C.check("analysis2_n", int(master["cc_analysis2"].sum()), label="Analysis 2 complete cases N")

    master.to_csv(C.MASTER, index=False)
    print(f"\n  -> {C.MASTER.name}  ({len(master)} rows x {master.shape[1]} columns)")


if __name__ == "__main__":
    sys.exit(main())
