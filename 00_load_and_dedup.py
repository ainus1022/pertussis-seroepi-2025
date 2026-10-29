"""
00_load_and_dedup.py - raw records -> curated, frozen participant files.

Outputs (private directory, participant-level):
  output/frozen_2025.csv     PWH 528 + PWoH 301 = 829
  output/frozen_2021.csv     PWH baseline (471 with anti-PT IgG)
  output/frozen_biobank.csv  serum-bank PWoH (149)

This step fixes the analysis sets. If the expected numbers cannot be reproduced
it stops. Dates of birth are used to verify identity but are never printed.
"""
import sys
import numpy as np
import pandas as pd

import config as C
import lib
from lib import to_num


# =============================================================================
# 2025 cohort (PWH + PWoH)
# =============================================================================
def load_2025() -> pd.DataFrame:
    d = pd.read_excel(C.SRC_2025_XLSX, sheet_name=C.SRC_2025_SHEET)
    d.columns = ["sample_no", "seq", "id_raw", "pt", "fha", "dob",
                 "age", "sex_ja", "dx", "cd4", "vl", "wbc", "lym"]

    d["pt"] = to_num(d["pt"])
    d["fha"] = to_num(d["fha"])
    d["seq"] = pd.to_numeric(d["seq"], errors="coerce")
    d["dob"] = pd.to_datetime(d["dob"])
    d["id"] = (d["id_raw"].astype(str).str.upper()
               .str.replace("-", "", regex=False).str.strip())

    n_raw = len(d)
    assert d["pt"].notna().all(), "2025 records: anti-PT IgG value that cannot be converted to a number"

    # --- Rule (3): correct a mistyped study ID (identified by date of birth) ---
    for wrong, dob, right in C.private("ID_TYPO_FIX"):
        m = (d["id"] == wrong) & (d["dob"] == pd.Timestamp(dob))
        if m.sum() != 1:
            raise AssertionError(f"ID correction: {m.sum()} matching records (expected exactly 1)")
        d.loc[m, "id"] = right
        print("  [rule 3] mistyped study ID corrected: 1 record")

    # --- Rule (4): documented transcription errors in a value (raw records are not edited) ---
    for vid, vdob, col, wrong_v, right_v in C.private("VALUE_FIX"):
        m = (d["id"] == vid) & (d[col] == wrong_v)
        if vdob is not None:
            m &= (d["dob"] == pd.Timestamp(vdob))
        if m.sum() != 1:
            raise AssertionError(f"value correction ({col}): {m.sum()} matching records (expected exactly 1)")
        d.loc[m, col] = right_v
        print(f"  [rule 4] transcription error in {col} corrected: 1 record")

    # --- Group membership comes from the study-ID prefix, not the diagnosis field ---
    d["group"] = np.where(d["id"].str.startswith(C.ID_PREFIX_PWH), "PWH", "PWoH")

    # --- Rule (2): one ID shared by two people -> keep both (different dates of birth) ---
    for i in C.private("KEEP_BOTH_IDS"):
        sub = d[d["id"] == i]
        if sub["dob"].nunique() != len(sub):
            raise AssertionError("shared ID: dates of birth are not distinct")
        print(f"  [rule 2] study ID shared by {len(sub)} different people: both kept")

    # --- Rule (1): one person measured twice -> keep the first measurement ---
    d = d.sort_values("seq").reset_index(drop=True)
    dropped = []
    for i in C.private("DUP_SAME_PERSON"):
        sub = d[d["id"] == i]
        if len(sub) != 2:
            raise AssertionError(f"repeat measurement: {len(sub)} records (expected 2)")
        if sub["dob"].nunique() != 1 or sub["sex_ja"].nunique() != 1:
            raise AssertionError("repeat measurement: date of birth or sex differs between records")
        dropped.extend(sub.index[1:])
        print("  [rule 1] participant measured twice: first measurement kept, second excluded")
    f = d.drop(index=dropped).reset_index(drop=True)

    # --- Final identity check: no ID + date-of-birth duplicate may remain ---
    key = f["id"] + "|" + f["dob"].dt.strftime("%Y%m%d")
    if key.duplicated().any():
        raise AssertionError(f"{int(key.duplicated().sum())} duplicate ID + date-of-birth records remain")

    f["male"] = (f["sex_ja"] == "男").astype(int)
    print(f"  2025: {n_raw} raw records -> {len(f)} after curation")
    return f


# =============================================================================
# Respiratory symptoms at blood draw (descriptive columns for the 2025 cohort)
# =============================================================================
def attach_symptoms(f: pd.DataFrame) -> pd.DataFrame:
    """Join the chart-review sheet to the frozen file by sample sequence number.

    The join key is the sequence number, not the study ID, because one study ID
    is shared by two people. Dates of birth must agree on every row afterwards.
    """
    sx = pd.read_excel(C.SRC_SYMPTOM_XLSX, sheet_name=C.SRC_SYMPTOM_SHEET)
    sx.columns = ["sample_no_sx", "draw_date", "sx_raw", "ims", "ka",
                  "dob_sx", "sex_sx", "dx_sx", "cd4_sx", "vl_sx", "wbc_sx", "lym_sx"]
    sx["seq"] = pd.to_numeric(
        sx["sample_no_sx"].astype(str).str.extract(r"D-(\d+)")[0], errors="coerce")
    if sx["seq"].isna().any() or sx["seq"].duplicated().any():
        raise AssertionError("symptom sheet: sequence number missing or duplicated")
    sx["dob_sx"] = pd.to_datetime(sx["dob_sx"])
    sx["sx_raw"] = sx["sx_raw"].astype(str).str.strip()

    m = f.merge(sx[["seq", "sample_no_sx", "draw_date", "sx_raw", "dob_sx"]],
                on="seq", how="left")
    if len(m) != len(f):
        raise AssertionError("symptom join changed the number of rows (key is not 1:1)")
    miss = int(m["sx_raw"].isna().sum())
    if miss:
        raise AssertionError(f"{miss} frozen records have no symptom record")
    bad = int((m["dob"] != m["dob_sx"]).sum())
    if bad:
        raise AssertionError(f"date of birth differs between the two sources on {bad} rows")

    m["any_sx"] = (m["sx_raw"] != C.SX_NONE).astype(int)
    m = m.drop(columns=["dob_sx"])
    print(f"  symptom sheet {len(sx)} rows -> joined to {len(m)} records (date of birth agrees on all)")
    print(f"  any symptom recorded: {int(m['any_sx'].sum())} ({100 * m['any_sx'].mean():.1f}%)"
          f"  [no note = not recorded, not absence of symptoms]")
    return m


# =============================================================================
# 2021 cohort (PWH, longitudinal baseline)
# =============================================================================
def load_2021() -> pd.DataFrame:
    d = pd.read_csv(C.SRC_2021_CSV)
    d = d.rename(columns={
        "サンプル no.": "sample_no", "生年月日": "dob_raw", "IMS no.": "id_raw",
        "PT-IgG": "pt", "FHA-IgG": "fha", "CD4数": "cd4", "CD8": "cd8",
        "HIV-RNA": "vl", "保存血清番号": "serum_no",
    })
    d["dob"] = pd.to_datetime(d["dob_raw"], format="mixed")
    d["id"] = (d["id_raw"].astype(str).str.upper()
               .str.replace("-", "", regex=False).str.strip())
    d["pt"] = to_num(d["pt"])

    n_raw = len(d)
    # Records without an anti-PT IgG result are outside the "PT cohort" but are kept:
    # the CD4/CD8 analysis uses the 2021 ratio as exposure and the 2025 titre as
    # outcome, so a 2021 titre is not required for it.
    f = d.copy()
    f["in_pt_cohort"] = f["pt"].notna().astype(int)

    # Identity check within the 2021 file.
    dup_id = f[f.duplicated("id", keep=False)]
    for _, sub in dup_id.groupby("id"):
        if sub["dob"].nunique() != len(sub):
            raise AssertionError("2021: duplicate ID with identical date of birth")
        print(f"  [rule 2] 2021: study ID shared by {len(sub)} different people: both kept")
    key = f["id"] + "|" + f["dob"].dt.strftime("%Y%m%d")
    if key.duplicated().any():
        raise AssertionError("2021: duplicate ID + date-of-birth records remain")

    f["cd48_ratio"] = f["cd4"] / f["cd8"]
    n_pt = int(f["in_pt_cohort"].sum())
    print(f"  2021: {n_raw} raw records -> {n_pt} with anti-PT IgG ({n_raw - n_pt} without)")
    print(f"        {int(f['cd48_ratio'].notna().sum())} with a CD4/CD8 ratio")
    return f


def main():
    print("=" * 72)
    print("00_load_and_dedup - raw records -> frozen participant files")
    print("=" * 72)

    print("\n[2025 cohort]")
    f25 = load_2025()
    n_pwh = int((f25["group"] == "PWH").sum())
    n_pwoh = int((f25["group"] == "PWoH").sum())
    C.check("n_pwh_2025", n_pwh, label="PWH 2025 n")
    C.check("n_pwoh_2025", n_pwoh, label="PWoH 2025 n")
    print(f"  OK  PWH={n_pwh} / PWoH={n_pwoh}")

    print("\n[Respiratory symptoms at blood draw]")
    f25 = attach_symptoms(f25)
    C.check("n_symptom_linked", len(f25), label="2025 records linked to a symptom record")
    C.check("n_any_sx", int(f25["any_sx"].sum()), label="any respiratory symptom recorded n")

    print("\n[2021 cohort]")
    f21 = load_2021()
    C.check("n_pwh_2021", int(f21["in_pt_cohort"].sum()), label="PWH 2021 n")
    print(f"  OK  PWH 2021={int(f21['in_pt_cohort'].sum())}")

    print("\n[Serum-bank PWoH]")
    bb = lib.load_biobank()
    C.check("n_biobank_pwoh_2025", len(bb), label="serum-bank PWoH n")
    print(f"  dispensed {C.BIOBANK['n_dispensed']} -> measured {len(bb)}")
    print(f"  OK  PWoH set A = {C.EXPECTED['n_pwoh_2025']} + {len(bb)} = "
          f"{C.EXPECTED['n_pooled_pwoh_2025']}; set B = {C.EXPECTED['n_pwoh_2025']}")

    cols25 = ["seq", "id", "group", "dob", "age", "male", "sex_ja", "dx",
              "pt", "fha", "cd4", "vl", "wbc", "lym",
              "sample_no_sx", "draw_date", "sx_raw", "any_sx"]
    f25[cols25].to_csv(C.FROZEN_2025, index=False)
    cols21 = ["sample_no", "id", "dob", "in_pt_cohort", "pt", "fha",
              "cd4", "cd8", "cd48_ratio", "vl", "serum_no"]
    f21[cols21].to_csv(C.FROZEN_2021, index=False)
    colsbb = ["code", "bb", "cohort", "birth_year", "age", "male", "sex",
              "age_grp_ja", "pt", "fha", "pt_pos", "pt_high"]
    bb[colsbb].to_csv(C.FROZEN_BIOBANK, index=False)
    print(f"\n  -> {C.FROZEN_2025.name}  ({len(f25)} rows)")
    print(f"  -> {C.FROZEN_2021.name}  ({len(f21)} rows)")
    print(f"  -> {C.FROZEN_BIOBANK.name}  ({len(bb)} rows)")


if __name__ == "__main__":
    sys.exit(main())
