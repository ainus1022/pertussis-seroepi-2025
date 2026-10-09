"""
lib.py - loaders shared by several scripts.

Numeric conversion and the construction of the PWoH analysis sets are defined
once here so that no script re-implements them.
"""
import numpy as np
import pandas as pd

import config as C


def to_num(s: pd.Series) -> pd.Series:
    """Convert assay values to numbers.

    Values at the upper limit of quantification are exported by the laboratory
    with a leading non-breaking space ('\\xa0160'). A plain pd.to_numeric turns
    them into missing values and undercounts seropositive participants, so
    non-breaking and ordinary spaces are stripped first.
    """
    return pd.to_numeric(
        s.astype(str).str.replace("\xa0", "", regex=False).str.strip(),
        errors="coerce",
    )


def vaccine_era(birth_year) -> str:
    """Birth year -> priming-vaccine generation (config section 11)."""
    if pd.isna(birth_year):
        return "unknown"
    y = int(birth_year)
    if y >= C.AP_BIRTH_YEAR_FROM:
        return "aP"
    if y >= C.WP_BIRTH_YEAR_UNTIL + 1:
        return "transition"
    return "wP"


def load_biobank() -> pd.DataFrame:
    """Load the serum-bank PWoH samples (healthy volunteers).

    A dispensed sample without an assay result is excluded (n = 149). Values
    above the upper limit of quantification are entered as 160 in the source
    workbook. Birth year is reconstructed as collection year minus age, because
    the serum bank provides sex and age only.
    """
    bb = pd.read_excel(C.SRC_BIOBANK_XLSX, sheet_name=C.SRC_BIOBANK_SHEET)
    bb = bb.rename(columns={"検体管理コード": "code", "BB番号": "bb", "性別": "sex",
                            "年齢": "age", "年代": "age_grp_ja", "採取日": "collected",
                            "PT-IgG(EU/mL)": "pt", "FHA-IgG(EU/mL)": "fha"})
    n_disp = len(bb)
    bb = bb[bb["pt"].notna()].copy()
    if n_disp != C.BIOBANK["n_dispensed"] or len(bb) != C.BIOBANK["n_analyzed"]:
        raise AssertionError(
            f"serum bank n: dispensed {n_disp} / measured {len(bb)} "
            f"(expected {C.BIOBANK['n_dispensed']}/{C.BIOBANK['n_analyzed']})")

    bb["pt"] = to_num(bb["pt"])
    bb["pt_pos"] = (bb["pt"] >= C.CUT_POS).astype(int)
    bb["pt_high"] = (bb["pt"] >= C.CUT_HIGH).astype(int)
    bb["birth_year"] = C.BIOBANK["collection_year"] - bb["age"]
    bb["male"] = (bb["sex"].astype(str).str.upper().str[0] == "M").astype(int)
    bb["cohort"] = "biobank"

    # Cross-check against the summary sheet of the source workbook.
    if int(bb["pt_pos"].sum()) != C.BIOBANK["pos_n"]:
        raise AssertionError(f"serum bank seropositive n {int(bb['pt_pos'].sum())} != {C.BIOBANK['pos_n']}")
    if int(bb["pt_high"].sum()) != C.BIOBANK["high_n"]:
        raise AssertionError(f"serum bank high-titre n {int(bb['pt_high'].sum())} != {C.BIOBANK['high_n']}")
    gmc = float(np.exp(np.log(bb["pt"].astype(float)).mean()))
    if abs(gmc - C.BIOBANK["gmc"]) > 0.15:
        raise AssertionError(f"serum bank GMC {gmc:.1f} != {C.BIOBANK['gmc']}")
    return bb


def _derive_pt_flags(df: pd.DataFrame) -> pd.DataFrame:
    """Derive the anti-PT IgG flags from the thresholds only.

    Derived columns in the raw or intermediate files are never trusted; flags are
    always recomputed from CUT_POS / CUT_HIGH / CUT_MULTI after pooling.
    """
    df = df.copy()
    df["pt_pos"] = (df["pt"] >= C.CUT_POS).astype(int)
    df["pt_high"] = (df["pt"] >= C.CUT_HIGH).astype(int)
    for cut in C.CUT_MULTI:
        df[f"pt_ge{cut}"] = (df["pt"] >= cut).astype(int)
    return df


def load_pwoh_set(caller: str, which: str | None = None) -> pd.DataFrame:
    """Return a PWoH analysis set (config section 10).

    A = hospital outpatients (301) + serum bank (149) = 450
    B = hospital outpatients only (301)

    The set is chosen by config.PWOH_SET_BY_ANALYSIS, not typed in by the caller.
    An explicit `which` is honoured, but a request for "A" always passes through
    the scope guard. Serum-bank rows have no wbc / lym / cd4 / vl / dx_cat; these
    are structurally missing and are not summarised for set A.
    """
    which = which or C.pwoh_set_for(caller)
    if which not in ("A", "B"):
        raise ValueError(f"load_pwoh_set: which='{which}' is invalid ('A' or 'B'). caller={caller}")

    m = pd.read_csv(C.MASTER, parse_dates=["dob"])
    main = m[m["hiv"] == 0].copy()
    main["birth_year"] = main["dob"].dt.year
    main["cohort"] = "main"

    if which == "B":
        out = _derive_pt_flags(main)
        if len(out) != C.EXPECTED["n_pwoh_2025"]:
            raise AssertionError(f"PWoH-B {len(out)} != {C.EXPECTED['n_pwoh_2025']}")
    else:
        C.require_pooled_scope(caller)
        bb = pd.read_csv(C.FROZEN_BIOBANK)
        # Union of columns; variables absent from the serum bank stay missing.
        out = _derive_pt_flags(pd.concat([main, bb], ignore_index=True))
        out["hiv"] = 0
        if len(out) != C.EXPECTED["n_pooled_pwoh_2025"]:
            raise AssertionError(f"PWoH-A {len(out)} != {C.EXPECTED['n_pooled_pwoh_2025']}")
        n_bb = int((out["cohort"] == "biobank").sum())
        if n_bb != C.EXPECTED["n_biobank_pwoh_2025"]:
            raise AssertionError(f"serum-bank component of PWoH-A {n_bb} != {C.EXPECTED['n_biobank_pwoh_2025']}")

    out["generation"] = out["birth_year"].map(vaccine_era)
    out.attrs["pwoh_set"] = which
    return out


def load_pooled_pwoh(caller: str) -> pd.DataFrame:
    """PWoH set A (n = 450); alias of load_pwoh_set(caller, "A")."""
    return load_pwoh_set(caller, "A")
