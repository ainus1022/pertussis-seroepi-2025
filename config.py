"""
config.py - single source of truth for the pertussis sero-epidemiology pipeline.

Every threshold, analysis-set rule, colour and expected value lives here. Each
script compares what it computes against EXPECTED and stops with an
ExpectationError on any mismatch, so a silent change in the data or the code
cannot propagate into the tables and figures.

Division of labour
  Python : cleaning, de-duplication, derived flags, every reported number and figure
  Stata  : independent recomputation of the test statistics only (verify_in_stata.do)

Privacy by construction
  Nothing that identifies a participant is stored in this repository. The raw
  input files, the participant-level intermediate files and the small set of
  record-level curation rules (which IDs are duplicates, hand-coded symptom
  categories) live in a separate private directory, by default ../private
  (override with the environment variable PERTUSSIS_PRIVATE_DIR). See
  config_private.example.py for the expected layout.
"""
import importlib.util
import os
from pathlib import Path

# =============================================================================
# 1. Paths
# =============================================================================

REPO = Path(__file__).resolve().parent
PRIVATE_DIR = Path(os.environ.get("PERTUSSIS_PRIVATE_DIR", REPO.parent / "private"))

DIR_NUM = REPO / "numbers"     # aggregate results reported in the paper (public)
DIR_FIG = REPO / "figures"     # figures (public)
DIR_DATA = REPO / "data"       # public reference data (JIHS 2023 survey table)
DIR_OUT = PRIVATE_DIR / "output"   # participant-level intermediates (never public)
for _d in (DIR_NUM, DIR_FIG, DIR_DATA):
    _d.mkdir(exist_ok=True)
if PRIVATE_DIR.exists():
    DIR_OUT.mkdir(exist_ok=True)


class PrivateDataUnavailable(RuntimeError):
    """Raised when a step needs the private directory and it is not present."""


def _load_private():
    f = PRIVATE_DIR / "config_private.py"
    if not f.exists():
        return None
    spec = importlib.util.spec_from_file_location("config_private", f)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_PRIV = _load_private()


def private(name: str):
    """Return a record-level curation rule from the private configuration."""
    if _PRIV is None or not hasattr(_PRIV, name):
        raise PrivateDataUnavailable(
            f"'{name}' is defined in the private configuration, which is not part of "
            f"this repository (expected at {PRIVATE_DIR / 'config_private.py'}). "
            f"Participant-level data are not publicly available; see README.md."
        )
    return getattr(_PRIV, name)


# --- Inputs (all private except the JIHS table). File names can be overridden
#     by an INPUTS dict in the private configuration.
INPUTS = {
    "cohort_2025": PRIVATE_DIR / "input" / "cohort_2025.xlsx",       # PWH + PWoH, 2025
    "cohort_2021": PRIVATE_DIR / "input" / "pwh_2021.csv",           # PWH baseline, 2021
    "serum_bank": PRIVATE_DIR / "input" / "serum_bank_2025.xlsx",    # serum-bank PWoH
    "symptoms": PRIVATE_DIR / "input" / "respiratory_symptoms.xlsx", # chart review
    "jihs_pdf": PRIVATE_DIR / "input" / "jihs_2023_survey.pdf",      # public report (optional)
    "cd8_2025": PRIVATE_DIR / "input" / "cd8_2025.csv",              # PWH CD8 at the 2025 sampling
}
if _PRIV is not None and hasattr(_PRIV, "INPUTS"):
    INPUTS.update({k: Path(v) for k, v in _PRIV.INPUTS.items()})

SRC_2025_XLSX = INPUTS["cohort_2025"]
SRC_2025_SHEET = "Sheet1"
SRC_2021_CSV = INPUTS["cohort_2021"]
SRC_BIOBANK_XLSX = INPUTS["serum_bank"]
SRC_BIOBANK_SHEET = "個票_n150"
SRC_SYMPTOM_XLSX = INPUTS["symptoms"]
SRC_SYMPTOM_SHEET = "呼吸器症状有無"
SRC_JIHS_PDF = INPUTS["jihs_pdf"]
SRC_CD8_2025_CSV = INPUTS["cd8_2025"]   # de-identified extract: id, dob, cd8 (other columns ignored)

# Age-specific anti-PT IgG distribution, general population, Japan, fiscal year 2023
# (National Epidemiological Surveillance of Vaccine-Preventable Diseases, Table 3-1).
# Extracted from the published report by 00b_extract_jihs.py and shipped with the
# repository because it is public aggregate data.
JIHS_FROM_PDF = DIR_DATA / "jihs2023_pt_from_pdf.csv"

FROZEN_2025 = DIR_OUT / "frozen_2025.csv"
FROZEN_2021 = DIR_OUT / "frozen_2021.csv"
FROZEN_BIOBANK = DIR_OUT / "frozen_biobank.csv"
MASTER = DIR_OUT / "analysis_master.csv"     # shared input for Python and Stata
POOLED_PWOH = DIR_OUT / "pooled_pwoh.csv"    # PWoH set A, input for the Stata check
KEY_NUMBERS = DIR_NUM / "key_numbers.json"
ENV_FILE = DIR_NUM / "env.txt"

# --- Figure file names. Change a name here and nowhere else.
FIGURES = {   # numbering follows the order of figures in the manuscript
    "flow":            "Fig1_participant_flow",    # participants and analysis sets
    "jihs_vs_pwoh":    "Fig2_jihs_vs_pwoh",        # general population 2023 vs PWoH 2025
    "age_pwoh_vs_pwh": "Fig3_age_pwoh_vs_pwh",     # seropositivity by age, PWoH vs PWH
    "forest":          "Fig4_forest",              # multivariable model
    "mcnemar_age":     "Fig5_mcnemar_age",         # PWH, paired 2021 vs 2025
    "cd48":            "Fig6_cd4cd8",              # PWH by baseline CD4/CD8 ratio
    "priming":         "Fig7_priming_generation",  # PWoH by priming-vaccine generation
}


def save_fig(fig, key: str, **kw):
    """Save a figure as PDF and PNG under the name registered in FIGURES."""
    if key not in FIGURES:
        raise KeyError(f"figure key not registered in FIGURES: {key}")
    stem = FIGURES[key]
    fig.savefig(DIR_FIG / f"{stem}.pdf", **kw)
    fig.savefig(DIR_FIG / f"{stem}.png", dpi=300, **kw)
    return stem


# =============================================================================
# 2. Assay and thresholds
# =============================================================================

# Anti-PT IgG and anti-FHA IgG were measured by a commercial laboratory
# (BML, Inc.; enzyme immunoassay) and are reported in EU/mL.
#
# EU/mL (assay-specific units) and IU/mL (international units) are not formally
# interchangeable. Values from other studies are therefore shown in the unit of
# the original publication and never converted (see TABLE2_UNIT_FOOTNOTE).
ASSAY_UNIT = "EU/mL"
LIT_UNIT = "IU/mL"
TABLE2_UNIT_FOOTNOTE = (
    "Titres from the present study are expressed in EU/mL as reported by the "
    "commercial laboratory (BML, Inc.; EIA), whereas the comparator studies report "
    "IU/mL. EU/mL and IU/mL are not formally interchangeable, and cut-offs are "
    "therefore not strictly comparable across rows; cross-study comparisons should "
    "be interpreted as approximate."
)
ASSAY_VENDOR = "BML, Inc."
ASSAY_METHOD = "EIA"

CUT_POS = 10                   # anti-PT IgG >= 10: seropositive (primary outcome)
CUT_HIGH = 100                 # anti-PT IgG >= 100: high titre (same definition as the JIHS survey)
CUT_MULTI = [10, 15, 20, 40]   # cut-offs used by the comparator studies in Table 2

# Upper limit of quantification (right-censoring). In the raw laboratory export
# these values carry a leading non-breaking space and must be cleaned before
# numeric conversion (lib.to_num), otherwise they are silently lost.
ASSAY_ULOQ = 160

# CD4/CD8 ratio strata: <= 0.5 versus > 0.5 (a ratio of exactly 0.5 is "low").
CD48_CUT = 0.5

# =============================================================================
# 3. Cohort definition
# =============================================================================
# Group membership is defined by the study-ID prefix (PWH versus PWoH series),
# not by the free-text diagnosis field.
#
# Three record-level curation rules are applied in 00_load_and_dedup.py. The
# rules are public; the specific IDs and dates of birth they apply to are in the
# private configuration (DUP_SAME_PERSON, KEEP_BOTH_IDS, ID_TYPO_FIX).
#   (1) A participant measured twice: keep the first measurement.
#   (2) One study ID shared by two different people: keep both
#       (confirmed as different people by date of birth).
#   (3) A mistyped study ID: corrected, confirmed by date of birth against the
#       2021 records. Without the correction one longitudinal pair is lost.
ID_PREFIX_PWH = "IMS"

# =============================================================================
# 4. Expected values
# =============================================================================
# "pwoh_*" without a suffix refers to PWoH set B (hospital outpatients);
# "pwoh_a_*" refers to set A (outpatients + serum bank). See section 10.

EXPECTED = {
    # --- Analysis sets (must never change) ---
    "n_pwh_2021": 471,
    "n_pwh_2025": 528,
    "n_pwoh_2025": 301,
    "n_biobank_pwoh_2025": 149,
    "n_pooled_pwoh_2025": 450,  # 301 + 149

    # --- Table 1 ---
    "pwh_age_median": 53,
    "pwoh_age_median": 63,
    "pwh_male_pct": 95.6,
    "pwh_lym_median": 1840, "pwoh_lym_median": 1610,
    "pwh_wbc_median": 5450,
    "pwh_cd4_median": 593.5,
    "pwh_vl_supp_n": 494, "pwh_vl_supp_denom": 524,   # HIV-RNA < 50 copies/mL
    "pwoh_liver_n": 139,
    "pwoh_metabolic_n": 54,
    "pwoh_checkup_n": 28,

    # --- Primary outcome ---
    "pwh_pos_n": 264, "pwh_pos_pct": 50.0,
    "pwoh_pos_n": 195, "pwoh_pos_pct": 64.8,
    "pwh_high_n": 13,
    "pwoh_high_n": 12,

    # --- PWoH set A ---
    "pwoh_a_n": 450,
    "pwoh_a_age_median": 55,
    "pwoh_a_male_n": 206, "pwoh_a_male_pct": 45.8,
    "pwoh_a_pos_n": 288, "pwoh_a_pos_pct": 64.0,
    "pwoh_a_high_n": 14,
    "pwoh_a_gmc": 15.5,
    "table1_a_pos_chi2": 19.368,
    "table1_a_age_p": 0.0268,

    # --- Analysis 1: general population (JIHS 2023, age >= 20) vs PWoH set A ---
    "jihs_pdf_n": 903, "jihs_pdf_pos_n": 449, "jihs_pdf_pos_pct": 49.7,
    "jihs_pdf_high_n": 4,
    "analysis1_a_pwoh_n": 450,
    "analysis1_a_pwoh_pos_n": 288,
    "analysis1_a_chi2_ge10": 24.684,
    "analysis1_a_chi2_ge50": 42.381,
    "analysis1_a_chi2_ge100": 16.288,

    # --- Analysis 2: multivariable logistic regression (outcome: anti-PT IgG >= 10) ---
    # Complete-case set = 829 participants in 2025 minus 96 with a missing
    # lymphocyte count = 733.
    "analysis2_n": 733,
    "analysis2_hiv_aor": 0.606,
    "analysis2_hiv_lci": 0.397,
    "analysis2_hiv_uci": 0.925,
    "analysis2_hiv_p": 0.020,

    # --- Analysis 2 sensitivity: age-matched pairs (MATCH_ALGORITHM) ---
    "matched_pairs": 163,
    "matched_or": 0.582,
    "matched_p": 0.018,

    # --- Analysis 3: PWH, paired 2021 -> 2025 ---
    "mcnemar_pairs": 463,
    "mcnemar_pos_2021": 225, "mcnemar_pos_2021_pct": 48.6,
    "mcnemar_pos_2025": 240, "mcnemar_pos_2025_pct": 51.8,
    "mcnemar_disc_neg_to_pos": 44,
    "mcnemar_disc_pos_to_neg": 29,
    # The reported statistic is McNemar's chi-square WITH continuity correction.
    # Stata's `mcc` reports the uncorrected statistic (3.082, p = 0.079); the two
    # packages agree on the exact test (p = 0.1006), which anchors the comparison.
    "mcnemar_chi2": 2.685,       # (|b-c|-1)^2/(b+c) = (15-1)^2/73
    "mcnemar_p": 0.101,
    "mcnemar_p_exact": 0.1006,
    "mcnemar_chi2_stata_uncorrected": 3.082,

    # --- Baseline (2021) CD4/CD8 ratio -> 2025 serostatus, PWH ---
    "cd48_pairs": 476,
    "cd48_low_n": 69, "cd48_low_pos": 27, "cd48_low_pos_pct": 39.1, "cd48_low_median": 7.0,
    "cd48_high_n": 407, "cd48_high_pos": 216, "cd48_high_pos_pct": 53.1, "cd48_high_median": 10.0,
    "cd48_fisher_p": 0.037,
    "cd48_mwu_p": 0.019,
    "cd48_exactly_half_n": 3,

    # --- Priming-vaccine generation, PWoH set A ---
    "priming_pooled_n": 450,
    "priming_ap_n": 125, "priming_ap_pos_pct": 54.4,
    "priming_wp_n": 273, "priming_wp_pos_pct": 68.5,
    "priming_transition_n": 52,
    "priming_model_n": 398,      # aP + wP (transition cohort excluded)
    "priming_or": 0.504,         # aP vs wP, adjusted for source cohort
    "priming_lci": 0.307, "priming_uci": 0.826, "priming_p": 0.0066,
    "priming_cohort_or": 1.201, "priming_cohort_p": 0.474,

    # --- Respiratory symptoms at blood draw (descriptive only) ---
    "n_symptom_linked": 829,
    "n_any_sx": 46,
    "n_sx_cough": 22,
    "n_sx_prolonged": 9,
    "n_sx_altdx": 7,

    # --- HIV disease history x 2025 seropositivity, within PWH (10_hiv_history.py) ---
    #     "CD4 at ART initiation" is the value recorded at ART start, not a true nadir.
    "hivhx_linked_n": 474,       # PWH 2025 linked to the 2021 baseline record sheet
    "hivhx_model_n": 443,        # complete cases for the multivariable model
    "hivhx_cd4now_p": 0.2158,    # chi2, <350 / 350-499 / >=500
    "hivhx_cd4art_p": 0.6955,    # chi2, <200 / 200-349 / >=350
    "hivhx_artdur_p": 0.1362,    # chi2, <10 / 10-19 / >=20 years
    "hivhx_aids_p": 0.8307,      # Fisher, yes / no
    # --- CD4/CD8 ratio at the 2025 sampling, within PWH (10_hiv_history.py) ---
    "cd48_2025_n": 521,          # PWH 2025 with a 2025 CD8 count (of 528)
    "cd48_2025_low_n": 74, "cd48_2025_low_pos_pct": 41.9,
    "cd48_2025_high_n": 447, "cd48_2025_high_pos_pct": 51.0,
    "cd48_2025_fisher_p": 0.168,
    "cd48_2025_mwu_p": 0.053,
    # --- priming-vaccine generation within PWH (aP born >=1982 vs wP born <=1974) ---
    "pwh_gen_ap_n": 106, "pwh_gen_ap_pos_pct": 40.6,
    "pwh_gen_wp_n": 317, "pwh_gen_wp_pos_pct": 52.1,
    "pwh_gen_trans_n": 105,
    "pwh_gen_fisher_p": 0.044,
    # --- raw row counts drawn in the participant flow (11_flowchart.py; figure input, not checked) ---
    "n_2025_raw": 831,           # 2025 records before de-duplication
    "n_2025_dup_removed": 2,     # second samples from the same person (831 -> 829)
    "n_2021_raw": 486,           # 2021 records, including 15 without anti-PT IgG (-> 471)
    "jihs_total_n": 1472,        # JIHS 2023 Table 3-1, all ages (adults: jihs_pdf_n)
}

# Tolerances: counts must match exactly; rounding differences are allowed elsewhere.
TOL_COUNT = 0
TOL_PCT = 0.15
TOL_STAT = 0.01

# =============================================================================
# 5. Respiratory symptoms at blood draw
# =============================================================================
# Source: retrospective review of outpatient charts (free text). Symptoms were not
# elicited systematically, so the absence of a note means "not recorded", not
# "no symptoms". The categories below were hand-coded AFTER the antibody results
# were known and are therefore used for description and for the Limitations only;
# they do not enter the primary or multivariable analyses.
#
# The hand-coded sets of sample numbers (SX_COUGH, SX_PROLONGED,
# SX_PROLONGED_STRICT, SX_ALTDX) are record-level and live in the private
# configuration.
#   SX_COUGH            cough or sputum recorded
#   SX_PROLONGED        cough recorded as lasting >= 2 weeks or as persistent
#   SX_PROLONGED_STRICT duration explicitly documented as >= 2 weeks
#   SX_ALTDX            another pathogen confirmed by testing
SX_NONE = "無"   # token used in the chart-review sheet for "no symptom recorded"

# =============================================================================
# 6. Analysis 2: covariates and the matched sensitivity analysis
# =============================================================================
# Logistic regression: outcome = anti-PT IgG >= CUT_POS; exposure = HIV;
# adjusted for age, sex and lymphocyte count. Complete cases only.
COVARIATES = ["age", "male", "lym"]
SCALE_AGE = 10       # odds ratio per 10 years
SCALE_LYM = 500      # odds ratio per 500 cells/uL

MATCH_AGE_CALIPER = 3
MATCH_ALGORITHM = (
    "1:1 nearest-neighbour matching on age within a ±3-year caliper, without replacement; "
    "PWH processed in ascending order of age (ties broken by participant ID); "
    "restricted to the Analysis-2 complete-case set; "
    "conditional odds ratio from discordant pairs (exact binomial test)."
)
# Sex is not a matching variable in the primary sensitivity analysis (it is adjusted
# for in the main model); a variant matched on age and sex is reported alongside.

# =============================================================================
# 7. Reason for hospital visit in PWoH (deterministic keyword classification)
# =============================================================================
# Priority: liver > metabolic > preventive > other. The keywords are matched
# against the Japanese diagnosis field of the hospital record.
DX_LIVER = [
    "肝炎", "肝機能障害", "脂肪肝", "NAFLD", "肝嚢胞", "肝腫瘍", "肝移植後",
    "薬剤性肝障害", "薬物性肝障害", "アルコール性肝障害",
    "原発性胆汁性胆管炎", "原発性胆汁性肝硬変", "原発性硬化性胆管炎",
    "ヘモクロマトーシス", "AFP高値",
]
DX_METABOLIC = [
    "高血圧", "糖尿病", "脂質異常症", "高脂血症", "高コレステロール血症", "高尿酸血症",
]
# Preventive visits: health check-up, vaccination, travel clinic, infection control.
DX_CHECKUP = [
    "健康診断", "健診", "渡航前相談", "渡航外来", "ワクチン接種", "院内感染対策",
]

# =============================================================================
# 8. Table 2: published comparators (cited by PMID)
# =============================================================================
TABLE2_PWOH = [
    # study, region (adult Tdap), n, age, cutoff, seropositive %, GMC, PMID
    ("Moriuchi 2017", "Japan (no adult Tdap)", "84", "4-44", "multilevel", "64.3", "58.4", "28700751"),
    ("Wehlin 2021", "Europe, 14 countries (Tdap, heterogeneous)", "7903", "20-39", "15", "53.6-88.0", "-", "34120372"),
    ("He 2022", "China (no)", "4459", "0-59", "20", "29.8", "9.67", "36283895"),
    ("Liu 2024 (Huzhou)", "China (no)", "1015", "0-83", "40", "11.43", "10.52", "38768133"),
    ("Greve-Isdahl 2025", "Norway (maternal programme; pregnant/mother-infant)", "366", "pregnant", "5 / 40", "52 / 5", "5.3", "41414933"),
    # The rows for the present study are generated from the data by 07_table2_international.py.
]

# High cut-offs reported for the present study in Table 2B:
# 94 as in Troy 2016, 100 as in Moosa 2025 and the JIHS survey.
TABLE2B_HIGH_CUTS = [94, 100]

TABLE2_PWH = [
    # study, region, n, age, CD4, cutoff, seropositive %, GMC, PMID
    ("Troy 2016", "USA", "299", "45.1 (mean)", "633 (mean)", "94", "4.3", "10.6", "27105050"),
    ("Moosa 2025", "South Africa (PHIRST; not HIV-specific)", "240 PWH / 1570", ">=65: 4.5%", "-", "100", "3.9-5.2", "-", "40645456"),
]

PMID = {
    "Moriuchi2017": "28700751",
    "Wehlin2021": "34120372",
    "He2022": "36283895",
    "Liu2024": "38768133",
    "Berbers2021": "34001895",
    "GreveIsdahl2025": "41414933",
    "Troy2016": "27105050",
    "Moosa2025": "40645456",
    "Dalby2010": "20508003",
    "deMelker2000": "10655388",
    "DeckerEdwards2021": "34590129",
    "Schwartz2016": "27672225",
    "SerranoVillar2014": "24831517",
    "Kimura1990": "1973217",     # introduction of acellular pertussis vaccine in Japan
    "Kanai1980": "7206322",      # 1975 suspension of vaccination and fall in coverage
    "Baughman2004": "15539504",  # source of the >= 94 diagnostic cut-off used by Troy 2016
}

# =============================================================================
# 9. Methods text
# =============================================================================
ETHICS = {
    "committee_en": "Research Ethics Committee, The Institute of Medical Science, The University of Tokyo",
    "approval_no": "2025-19-0722",
    "consent": "opt-out",
}

# The Stata step recomputes test statistics from the same analysis file; it does
# not re-derive the cohort. The wording is deliberately limited to that.
METHODS_STATS_SENTENCE = (
    "All analyses were performed using Python {py} with pandas {pd}, NumPy {np}, "
    "SciPy {sp}, and statsmodels {sm}; figures were generated with Matplotlib {mpl}. "
    "The reported test statistics were independently recomputed in Stata (version 19; "
    "StataCorp, College Station, TX, USA)."
)


class ExpectationError(AssertionError):
    """A computed value differs from EXPECTED. run_all.py stops here."""


_LOG: list[tuple[str, object, object, bool]] = []


def check(key: str, actual, tol: float | None = None, label: str | None = None):
    """Compare `actual` with EXPECTED[key]; raise ExpectationError on mismatch."""
    if key not in EXPECTED:
        raise KeyError(f"key not defined in EXPECTED: {key}")
    exp = EXPECTED[key]
    if tol is None:
        tol = TOL_COUNT if isinstance(exp, int) else TOL_STAT
    ok = abs(float(actual) - float(exp)) <= tol
    _LOG.append((label or key, exp, actual, ok))
    if not ok:
        raise ExpectationError(
            f"[{label or key}] expected={exp} actual={actual} (tol={tol})\n"
            f"  Either update config.EXPECTED or fix the derivation - decide explicitly."
        )
    return actual


def check_pct(key: str, actual, label: str | None = None):
    return check(key, actual, tol=TOL_PCT, label=label)


def checks_log():
    return list(_LOG)


# =============================================================================
# 10. PWoH analysis sets A and B
# =============================================================================
# A = hospital outpatients (301) + IMSUT serum bank (149) = 450
#     The main PWoH set. The serum-bank samples come from healthy volunteers
#     (additional blood draw at staff health check-ups) and carry sex and age only.
# B = hospital outpatients only = 301
#     The subset with blood counts and reason for visit.
#
# The multivariable model adjusts for lymphocyte count, which does not exist for
# serum-bank samples, so it is necessarily restricted to B. Every other PWoH
# analysis needs only anti-PT IgG, age and sex and uses A. Table 1 shows both.

BIOBANK = {
    "n_dispensed": 150,
    "n_analyzed": 149,          # one dispensed sample has no assay result
    "collection_year": 2025,
    # Totals from the summary sheet of the source workbook (transcription check).
    "pos_n": 93, "pos_pct": 62.4, "high_n": 2, "gmc": 14.7,
    "source_name_en": "IMSUT serum bank",
    "source_url": "https://www.ims.u-tokyo.ac.jp/serumbank/index.html",
    "population": "healthy volunteers (additional blood draw at staff health check-ups)",
    # Sex and age only; birth year is reconstructed as collection year minus age
    # (error of up to one year, absorbed by excluding the transition cohort).
    "available_covariates": ["sex", "age"],
    "sampling": "all available samples provided (convenience)",
}

BIOBANK_ACKNOWLEDGMENT_EN = (
    "Serum samples from adults without HIV were provided by the IMSUT serum bank, "
    "The Institute of Medical Science, The University of Tokyo."
)

PWOH_SETS = {
    "A": "IMSUT Hospital outpatients + IMSUT serum bank (pooled)",   # 450
    "B": "IMSUT Hospital outpatients only",                          # 301
}

# Which PWoH set each script uses: "A" / "B" / "AB" (both shown) / "none" / "biobank".
PWOH_SET_BY_ANALYSIS = {
    "02_table1_characteristics":  "AB",
    "03_analysis1_jihs_vs_pwoh":  "A",
    "04_analysis2_forest":        "B",        # needs lymphocyte count
    "05_analysis3_mcnemar":       "none",     # PWH only
    "06_cd4cd8_stratified":       "none",     # PWH only
    "07_table2_international":    "A",
    "08_biobank_descriptive":     "biobank",
    "09_fig6_priming_generation": "A",
}

# Variables that cannot be reported for set A (structurally missing in the serum bank).
PWOH_A_UNAVAILABLE = ["wbc", "lym", "cd4", "vl", "dx_cat"]
PWOH_A_UNAVAILABLE_NOTE = (
    "Not available for dataset A: the serum-bank component provides sex and age only, "
    "with no haematology or clinical indication. These variables are shown for dataset B "
    "(hospital outpatients), which is the set used for the multivariable model."
)


def pwoh_set_for(caller: str) -> str:
    """Return the PWoH set assigned to a script (stop if it is not registered)."""
    if caller not in PWOH_SET_BY_ANALYSIS:
        raise KeyError(
            f"no PWoH set assigned to {caller}; add it to config.PWOH_SET_BY_ANALYSIS."
        )
    return PWOH_SET_BY_ANALYSIS[caller]


def require_pooled_scope(caller: str):
    """Guard against a set-B analysis accidentally loading the pooled set A."""
    which = pwoh_set_for(caller)
    if which in ("A", "AB"):
        return
    raise PermissionError(
        f"{caller} is assigned PWoH set '{which}' and must not load the pooled set A "
        f"(n={EXPECTED['n_pooled_pwoh_2025']})."
    )


# =============================================================================
# 11. Priming-vaccine generation (by birth year)
# =============================================================================
#   aP         : born >= 1982       acellular vaccine era
#   transition : born 1975-1981     low and heterogeneous coverage after the 1975 suspension
#   wP         : born <= 1974       whole-cell vaccine era
# Sources: Kimura & Kuno-Sakai, Lancet 1990 (PMID 1973217); Kanai, Jpn J Med Sci
# Biol 1980 (PMID 7206322); JIHS survey report (DTaP in routine use since 1982).
AP_BIRTH_YEAR_FROM = 1982
WP_BIRTH_YEAR_UNTIL = 1974

# Sensitivity 1: two groups including the transition cohort, split at 1982 or 1981.
PRIMING_SENSITIVITY_CUTS = [1982, 1981]
# Sensitivity 2: lower bound of the transition band.
PRIMING_TRANSITION_LOWER = [1975, 1974, 1973]

# =============================================================================
# 12. Colours: one colour per population in every figure (Okabe-Ito based)
# =============================================================================
PALETTE = {
    "pwoh_2025": "#0072B2",   # PWoH 2025
    "jihs_2023": "#8C8C8C",   # general population 2023 (reference, neutral grey)
    "pwh_2025":  "#D55E00",   # PWH 2025
    "pwh_2021":  "#6A3D9A",   # PWH 2021
    "wp":        "#009E73",   # wP generation
    "transition":"#56B4E9",   # transition cohort
    "ap":        "#CC79A7",   # aP generation
    "ink":       "#22303F",   # point estimates, median lines
    "ref_line":  "#9AA5B1",   # reference lines (OR = 1, cut-off)
}
