# Pertussis toxin IgG seropositivity during the 2025 outbreak in Japan: analysis code

Analysis code and aggregate results for

> Okazaki A, Koga M, et al. **Pertussis Toxin IgG Seropositivity During the 2025
> Outbreak: A Comparison Between People With HIV and Without HIV.**
> IDWeek 2026, Washington, DC. Poster P-1187.

A single-centre sero-epidemiological study at IMSUT Hospital, The Institute of
Medical Science, The University of Tokyo. Anti-pertussis toxin (PT) IgG was
measured in residual sera from people with HIV (PWH) and adults without HIV
(PWoH) during the 2025 pertussis epidemic, and in paired sera of PWH from 2021.

The same code will accompany the full manuscript; a tagged release and DOI will
be added here at submission.

## What is in this repository

| Path | Content |
|---|---|
| `config.py` | Single source of truth: thresholds, analysis-set rules, expected values, colours |
| `00_`–`09_*.py`, `lib.py`, `run_all.py` | The analysis pipeline |
| `poster_figures.py` | Figures as laid out on the IDWeek 2026 poster |
| `symptoms_descriptive.py` | Description of respiratory symptoms at blood draw |
| `verify_in_stata.do` | Independent recomputation of the test statistics in Stata |
| `numbers/` | Every reported number, as aggregate CSV files |
| `figures/` | Manuscript figures; `figures/poster/` poster figures |
| `data/jihs2023_pt_from_pdf.csv` | General-population reference (public aggregate data, see below) |

## What is not in this repository

Participant-level data are not publicly available to protect participant
privacy. The study uses residual clinical sera under an opt-out consent
procedure approved by the Research Ethics Committee of The Institute of Medical
Science, The University of Tokyo (approval no. 2025-19-0722).

The pipeline is built so that nothing identifying is ever written inside this
repository. Raw inputs, participant-level intermediate files, and the few
record-level curation rules (which study IDs are repeat measurements, the
hand-coded symptom categories) live in a separate private directory
(`../private` by default, or `$PERTUSSIS_PRIVATE_DIR`). The layout is described
in `config_private.example.py`. Without that directory the code can be read and
the aggregate results inspected, but the pipeline cannot be re-run.

## Pipeline

| Script | Step | Output |
|---|---|---|
| `00_load_and_dedup.py` | Raw records → curated participant files; identity checks by date of birth | private |
| `00b_extract_jihs.py` | Extracts or validates the general-population reference table | `data/` |
| `01_derive_and_validate.py` | Derives every flag from the thresholds; builds the single analysis file | private |
| `02_table1_characteristics.py` | Table 1; seropositivity by age group, PWoH vs PWH | `table1*.csv`, Fig2 |
| `03_analysis1_jihs_vs_pwoh.py` | General population 2023 vs PWoH 2025 | `analysis1*.csv`, Fig1 |
| `04_analysis2_forest.py` | Multivariable logistic regression; age-matched sensitivity analysis | `analysis2*.csv`, Fig3 |
| `05_analysis3_mcnemar.py` | PWH, paired sera 2021 → 2025 | `mcnemar_*.csv`, Fig4 |
| `06_cd4cd8_stratified.py` | PWH by baseline CD4/CD8 ratio | `cd4cd8.csv`, Fig5 |
| `07_table2_international.py` | Comparison with published studies (cited by PMID) | `table2*.csv` |
| `08_biobank_descriptive.py` | Serum-bank component of the PWoH set | `biobank_n149.csv` |
| `09_fig6_priming_generation.py` | PWoH by priming-vaccine generation (whole-cell vs acellular) | `priming_*.csv`, Fig6 |
| `10_hiv_history.py` | PWH: CD4 at ART initiation, ART duration, AIDS history and current CD4 vs 2025 seropositivity | `hiv_history*.csv` |

Figure file names carry the figure numbers of the manuscript (Fig1 general
population vs PWoH; Fig2 PWoH vs PWH by age; Fig3 multivariable model; Fig4
paired 2021 vs 2025; Fig5 baseline CD4/CD8 ratio; Fig6 priming-vaccine generation).

Poster figures (`poster_figures.py`): `Poster_Fig1A` general population vs PWoH
by age; `1B` PWoH by priming-vaccine generation; `2A` PWoH vs PWH by age; `2B`
multivariable model; `3A` PWH paired 2021 vs 2025; `3B` 2025 titre by baseline
CD4/CD8 ratio.

## Analysis sets

| Set | n | Definition |
|---|---|---|
| PWH 2025 | 528 | People with HIV attending IMSUT Hospital, 2025 |
| PWH 2021 | 471 | Stored sera from 2021 with an anti-PT IgG result |
| PWoH set A | 450 | Hospital outpatients without HIV (301) + IMSUT serum bank, healthy volunteers (149) |
| PWoH set B | 301 | Hospital outpatients without HIV only |
| General population 2023 | 903 | JIHS survey, age ≥ 20 |
| Multivariable model | 733 | PWH + PWoH set B with age, sex and lymphocyte count available |
| Paired PWH | 463 | Anti-PT IgG in both 2021 and 2025 |
| CD4/CD8 analysis | 476 | 2021 CD4/CD8 ratio and 2025 anti-PT IgG |

Set A is the main PWoH set. Serum-bank samples come with sex and age only, so
the multivariable model, which adjusts for lymphocyte count, is necessarily
restricted to set B. Table 1 reports both.

Seropositivity is anti-PT IgG ≥ 10 EU/mL and high titre is ≥ 100 EU/mL (enzyme
immunoassay, BML, Inc.). EU/mL are assay-specific units and are not converted to
IU/mL; values from other studies are shown in their published units.

## Reproducibility

Every script compares its results with `config.EXPECTED` and stops on any
difference, so the numbers in `numbers/` cannot drift from the code silently.

```bash
pip install -r requirements.txt        # pdftotext (poppler) is needed only to re-extract the JIHS table
python run_all.py --biobank --symptoms
python poster_figures.py
```

`numbers/key_numbers.json` records the main results and the software versions
of the run that produced them.

Test statistics were independently recomputed in Stata 19 from the same analysis
file (`verify_in_stata.do`). This checks the tests, not the derivation of the
cohort. McNemar's test is reported with continuity correction; Stata's `mcc`
reports the uncorrected statistic, and the two programs agree on the exact test.

## General-population reference data

`data/jihs2023_pt_from_pdf.csv` reproduces Table 3-1 (age-specific distribution
of anti-PT IgG) of the fiscal-year 2023 report of the National Epidemiological
Surveillance of Vaccine-Preventable Diseases, published by the Japan Institute
for Health Security (JIHS). It is extracted from the published report by
`00b_extract_jihs.py`, which verifies row and column totals against the report.

## Acknowledgement

Serum samples from adults without HIV were provided in part by the IMSUT serum
bank, The Institute of Medical Science, The University of Tokyo.

## Licence and contact

Code: MIT licence (see `LICENSE`). Aggregate results and figures: CC BY 4.0.

Aiko Okazaki, MD, PhD — UTOPIA and The Institute of Medical Science, The
University of Tokyo — okazaki-aiko321@g.ecc.u-tokyo.ac.jp
