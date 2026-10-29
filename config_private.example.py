"""
config_private.example.py - template for the private configuration.

The real file is named config_private.py and lives OUTSIDE this repository, in
the private directory (default ../private, or $PERTUSSIS_PRIVATE_DIR):

    private/
        config_private.py      record-level curation rules (this template)
        input/                 raw input files
        output/                participant-level intermediate files (written by the pipeline)

All values below are placeholders. The real values identify individual records
and are not publicly available.
"""

# Optional: locations of the raw input files, if they are not in private/input/.
# INPUTS = {
#     "cohort_2025": "/path/to/cohort_2025.xlsx",
#     "cohort_2021": "/path/to/pwh_2021.csv",
#     "serum_bank":  "/path/to/serum_bank_2025.xlsx",
#     "symptoms":    "/path/to/respiratory_symptoms.xlsx",
#     "jihs_pdf":    "/path/to/jihs_2023_survey.pdf",
#     "cd8_2025":    "/path/to/cd8_2025.csv",   # columns: id, dob, cd8 (de-identified extract)
# }

# Rule (1): study IDs of participants measured twice (first measurement kept).
DUP_SAME_PERSON = ["XX000", "YYY0000"]

# Rule (2): study IDs shared by two different people (both kept).
KEEP_BOTH_IDS = ["YYY0001"]

# Rule (3): mistyped study IDs: (ID as recorded, date of birth, correct ID).
ID_TYPO_FIX = [
    ("YYY0002", "1900-01-01", "YYY0003"),
]

# Rule (4): documented transcription errors in a value:
#   (ID, date of birth or None, column, wrong value, corrected value).
VALUE_FIX = [
    ("YYY0004", None, "cd4", 1000, 100),
]

# Hand-coded respiratory-symptom categories (sample numbers from the chart review).
SX_COUGH = {"D-000", "D-001"}             # cough or sputum recorded
SX_PROLONGED = {"D-000"}                  # cough >= 2 weeks or persistent (subset of SX_COUGH)
SX_PROLONGED_STRICT = {"D-000"}           # duration explicitly >= 2 weeks (subset of SX_PROLONGED)
SX_ALTDX = {"D-002"}                      # another pathogen confirmed by testing

assert SX_PROLONGED <= SX_COUGH
assert SX_PROLONGED_STRICT <= SX_PROLONGED
