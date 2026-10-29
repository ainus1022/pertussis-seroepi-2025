*==============================================================================
* verify_in_stata.do - independent recomputation of the test statistics
*                      (Stata 19, StataCorp)
*
* What this checks   : the test statistics themselves (aOR, chi-square, p),
*                      i.e. that the Python calls use the intended test
*                      (two-sided, tie handling, continuity correction).
* What it does not   : the derivation of the cohort (de-duplication, thresholds,
*                      flags). Both programs read the same analysis file.
*
* Usage: set the working directory to the folder containing this file, then
*          do verify_in_stata.do
*        The private directory defaults to ../private; to use another location:
*          do verify_in_stata.do "/path/to/private"
*==============================================================================

args priv
if "`priv'" == "" local priv "../private"

capture log close
log using "`priv'/output/verify_stata_log.txt", replace text

version 19
clear all
set more off

import delimited "`priv'/output/analysis_master.csv", clear varnames(1) case(preserve)

di as txt "{hline 78}"
di as res "0. Analysis sets (expected PWH 528 / PWoH 301)"
di as txt "{hline 78}"
tab group

*------------------------------------------------------------------------------
* 1. Table 1, PWH vs PWoH set B (hospital outpatients)
*    Expected: seropositive 50.0% (264/528) vs 64.8% (195/301), p<0.001
*------------------------------------------------------------------------------
di as res "1. Table 1: PWH vs PWoH set B"
tab hiv pt_pos, chi2 row
ranksum age, by(hiv)
ranksum lym, by(hiv)
ranksum wbc, by(hiv)

*------------------------------------------------------------------------------
* 2. Analysis 2: multivariable logistic regression (complete cases)
*    Expected: N=733, HIV aOR 0.606 (0.397-0.925) p=0.020
*              age (per 10 y) 1.048 / male 0.967 / lymphocytes (per 500) 0.958
*------------------------------------------------------------------------------
di as txt "{hline 78}"
di as res "2. Analysis 2: multivariable logistic regression"
di as txt "{hline 78}"
gen age_s = age/10
gen lym_s = lym/500
logit pt_pos hiv age_s male lym_s if cc_analysis2==1
logistic pt_pos hiv age_s male lym_s if cc_analysis2==1
estat ic

*------------------------------------------------------------------------------
* 3. Analysis 3: McNemar (PWH, paired 2021 -> 2025)
*    Expected: 463 pairs, 48.6% (225) -> 51.8% (240); discordant 44 vs 29
*    Stata's mcc/mcci report the UNCORRECTED statistic:
*        (b-c)^2/(b+c)     = (44-29)^2/73 = 3.08,  p = 0.0792
*    The reported value uses the continuity correction (Python, statsmodels):
*        (|b-c|-1)^2/(b+c) = (15-1)^2/73  = 2.685, p = 0.1013
*    The exact test agrees between the two programs (p = 0.1006).
*------------------------------------------------------------------------------
di as txt "{hline 78}"
di as res "3. Analysis 3: McNemar"
di as txt "{hline 78}"
mcc pt_pos_2021 pt_pos if paired==1
mcci 196 29 44 194

*------------------------------------------------------------------------------
* 4. Baseline (2021) CD4/CD8 ratio -> 2025 serostatus
*    Expected: n=476; <=0.5 N=69, 39.1%; >0.5 N=407, 53.1%
*              Fisher p=0.037 / Mann-Whitney p=0.019
*    Stata's ranksum applies no continuity correction; SciPy's mannwhitneyu does
*    by default. Compare p-values, not the U statistic.
*------------------------------------------------------------------------------
di as txt "{hline 78}"
di as res "4. CD4/CD8 ratio strata"
di as txt "{hline 78}"
gen byte cd48_low = (cd48_ratio_2021 <= 0.5) if cd48_linked==1
tab cd48_low pt_pos if cd48_linked==1, exact chi2 row
ranksum pt if cd48_linked==1, by(cd48_low)

*------------------------------------------------------------------------------
* 5. PWoH set A (outpatients + serum bank, n=450)
*    5a. Table 1: PWH vs PWoH set A
*        Expected: 264/528 (50.0%) vs 288/450 (64.0%), chi2 = 19.368
*    5b. Analysis 1: general population (JIHS 2023, age >= 20) vs PWoH set A
*        Expected: >=10   449/903 vs 288/450, chi2 = 24.684
*                  >=50    33/903 vs  59/450, chi2 = 42.381
*                  >=100    4/903 vs  14/450, chi2 = 16.288
*        The survey provides counts only, so each 2x2 table is entered as counts
*        and expanded to one row per person.
*    5c. Priming-vaccine generation
*        Expected: wP n=273 68.5% / transition n=52 63.5% / aP n=125 54.4%
*                  aP vs wP (N=398, adjusted for source cohort):
*                    aOR 0.504 (0.307-0.826) p=0.0066; cohort OR 1.201 p=0.474
*                  crude OR 0.549, Fisher p=0.0071
*------------------------------------------------------------------------------
di as txt "{hline 78}"
di as res "5a. Table 1: PWH vs PWoH set A"
di as txt "{hline 78}"
preserve
clear
input byte pwoh byte pos int freq
0 1 264
0 0 264
1 1 288
1 0 162
end
expand freq
label define grpA 0 "PWH 2025" 1 "PWoH-A 2025"
label values pwoh grpA
tab pwoh pos, chi2 exact row
restore

di as txt "{hline 78}"
di as res "5b. Analysis 1: JIHS 2023 vs PWoH set A"
di as txt "{hline 78}"
preserve

clear
input byte pwoh byte pos int freq
0 1 449
0 0 454
1 1 288
1 0 162
end
expand freq
label define grp 0 "JIHS 2023" 1 "PWoH-A 2025"
label values pwoh grp
di as res "  [>=10]"
tab pwoh pos, chi2 exact row

clear
input byte pwoh byte mid int freq
0 1 33
0 0 870
1 1 59
1 0 391
end
expand freq
label values pwoh grp
di as res "  [>=50]"
tab pwoh mid, chi2 exact row

clear
input byte pwoh byte high int freq
0 1 4
0 0 899
1 1 14
1 0 436
end
expand freq
label values pwoh grp
di as res "  [>=100]"
tab pwoh high, chi2 exact row

restore

di as txt "{hline 78}"
di as res "5c. Priming-vaccine generation (PWoH set A)"
di as txt "{hline 78}"
preserve
import delimited "`priv'/output/pooled_pwoh.csv", clear varnames(1) case(preserve)
tab cohort pt_pos, row
tab generation cohort
tab generation pt_pos, row

gen byte aP = (generation == "aP") if inlist(generation, "aP", "wP")
capture drop bb   // the file already carries bb; rebuild it from cohort
gen byte bb = (cohort == "biobank")
tab aP pt_pos if !missing(aP), exact chi2 row
logistic pt_pos aP bb if !missing(aP)

* Sensitivity: two groups including the transition cohort (split at 1982 / 1981)
gen byte aP82 = (birth_year >= 1982)
gen byte aP81 = (birth_year >= 1981)
tab aP82 pt_pos, exact row
tab aP81 pt_pos, exact row
restore

di as txt "{hline 78}"
di as res "6. HIV disease history -> 2025 serostatus (PWH only; 10_hiv_history.py)"
di as txt "{hline 78}"
* Expected (numbers/hiv_history.csv, hiv_history_model.csv):
*   CD4 2025       <350 24/58 41.4% | 350-499 50/109 45.9% | >=500 186/357 52.1% | chi2 p=0.216
*   CD4 at ART     <200 105/203 51.7% | 200-349 92/173 53.2% | >=350 34/72 47.2%  | chi2 p=0.696
*   ART duration   <10 31/74 41.9% | 10-19 149/289 51.6% | >=20 62/109 56.9%      | chi2 p=0.136
*   AIDS history   no 180/352 51.1% | yes 61/116 52.6%                            | Fisher p=0.831
*   logistic (n=443): age/10 1.015, male 1.006, cd4_2025/100 1.044, cd4_art/100 0.970,
*                     art/5y 1.110 (0.931-1.323, p=0.245), aids 1.141
preserve
import delimited "`priv'/output/hiv_history_cc.csv", clear varnames(1) case(preserve)
tab cd4_now_cat pt_pos, chi2 row
tab cd4_art_cat pt_pos, chi2 row
tab art_cat pt_pos, chi2 row
tab aids pt_pos, exact row
gen double age10 = age/10
gen double cd4_100 = cd4/100
gen double cd4art_100 = cd4_art/100
gen double art5 = art_years/5
logistic pt_pos age10 male cd4_100 cd4art_100 art5 aids
restore

* 6b. 2025 CD4/CD8 ratio and priming generation within PWH
*     Expected: ratio<=0.5 31/74 41.9% vs >0.5 228/447 51.0%, Fisher p=0.168, ranksum p~0.053
*               logistic per unit ratio (age, sex adj.): OR 1.751 (1.195-2.566) p=0.004
*               aP 43/106 40.6% vs wP 165/317 52.1%, Fisher p=0.044
preserve
import delimited "`priv'/output/hiv_history_cc.csv", clear varnames(1) case(preserve)
gen byte low25 = (cd48_ratio_2025 <= 0.5) if !missing(cd48_ratio_2025)
tab low25 pt_pos, exact row
ranksum pt, by(low25)
gen double age10 = age/10
logistic pt_pos age10 male cd48_ratio_2025
tab generation pt_pos, row
tab generation pt_pos if generation != "transition", exact row
restore

di as txt "{hline 78}"
di as res "Done. Compare with numbers/*.csv."
di as txt "{hline 78}"

log close
