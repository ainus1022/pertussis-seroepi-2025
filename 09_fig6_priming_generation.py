"""
09_fig6_priming_generation.py - PWoH: seropositivity by priming-vaccine generation.

Analysis set: PWoH set A = hospital outpatients (301) + serum bank (149) = 450.
The outpatients alone are older and contain few adults primed with acellular
vaccine; the serum bank supplies the younger age range.

Main analysis : aP (born >= 1982) vs wP (born <= 1974); transition cohort
                (1975-1981) excluded; adjusted for source cohort.
Sensitivity 1 : two groups including the transition cohort (split at 1982 / 1981).
Sensitivity 2 : lower bound of the transition band varied.

Input : output/analysis_master.csv, output/frozen_biobank.csv   (private)
Output: numbers/priming_generation.csv, numbers/priming_sensitivity.csv,
        numbers/priming_transition_band.csv, figures/Fig6_*
"""
import sys
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import config as C
import lib

ORDER = ["wP", "transition", "aP"]
LABEL = {
    "wP": f"wP\n(born $\\leq${C.WP_BIRTH_YEAR_UNTIL})",
    "transition": f"Transition\n({C.WP_BIRTH_YEAR_UNTIL + 1}-{C.AP_BIRTH_YEAR_FROM - 1})",
    "aP": f"aP\n(born $\\geq${C.AP_BIRTH_YEAR_FROM})",
}


def gm(s):
    return float(np.exp(np.log(s.astype(float)).mean()))


def main():
    print("=" * 72)
    print("09_fig6_priming_generation - priming-vaccine generation (PWoH)")
    print("=" * 72)

    p = lib.load_pooled_pwoh("09_fig6_priming_generation")
    C.check("priming_pooled_n", len(p), label="PWoH set A n")
    # Input for the Stata check (the analysis file does not contain the serum bank).
    p.to_csv(C.POOLED_PWOH, index=False)
    print(f"  PWoH set A n={len(p)} "
          f"(outpatients {int((p['cohort'] == 'main').sum())} + serum bank {int((p['cohort'] == 'biobank').sum())})")

    # --- Generation by source cohort (shows why adjustment for cohort is needed) ---
    cross = pd.crosstab(p["generation"], p["cohort"]).reindex(ORDER)
    print("\n[Generation x source cohort]")
    print(cross.to_string())

    # --- Summary by generation ---
    rows = []
    for g in ORDER:
        s = p[p["generation"] == g]
        rows.append({
            "generation": g,
            "birth_years": (f"<={C.WP_BIRTH_YEAR_UNTIL}" if g == "wP" else
                            f">={C.AP_BIRTH_YEAR_FROM}" if g == "aP" else
                            f"{C.WP_BIRTH_YEAR_UNTIL + 1}-{C.AP_BIRTH_YEAR_FROM - 1}"),
            "N": len(s), "seropositive_n": int(s["pt_pos"].sum()),
            "seropositive_pct": round(100 * s["pt_pos"].mean(), 1),
            "GMC": round(gm(s["pt"]), 1),
            "age_median": round(s["age"].median(), 0),
            "n_main": int((s["cohort"] == "main").sum()),
            "n_biobank": int((s["cohort"] == "biobank").sum()),
        })
    summ = pd.DataFrame(rows)
    C.check("priming_ap_n", int(summ.loc[summ.generation == "aP", "N"].iloc[0]), label="aP n")
    C.check("priming_wp_n", int(summ.loc[summ.generation == "wP", "N"].iloc[0]), label="wP n")
    C.check("priming_transition_n",
            int(summ.loc[summ.generation == "transition", "N"].iloc[0]), label="transition n")
    C.check_pct("priming_ap_pos_pct",
                float(summ.loc[summ.generation == "aP", "seropositive_pct"].iloc[0]), label="aP seropositive %")
    C.check_pct("priming_wp_pos_pct",
                float(summ.loc[summ.generation == "wP", "seropositive_pct"].iloc[0]), label="wP seropositive %")

    # --- Main analysis: aP vs wP, transition excluded, adjusted for source cohort ---
    sub = p[p["generation"].isin(["aP", "wP"])].copy()
    sub["aP"] = (sub["generation"] == "aP").astype(int)
    sub["biobank"] = (sub["cohort"] == "biobank").astype(int)
    C.check("priming_model_n", len(sub), label="priming model N")
    r = sm.Logit(sub["pt_pos"], sm.add_constant(sub[["aP", "biobank"]].astype(float))).fit(disp=0)
    ci = r.conf_int()
    or_ap = float(np.exp(r.params["aP"]))
    lci, uci = float(np.exp(ci.loc["aP", 0])), float(np.exp(ci.loc["aP", 1]))
    or_co = float(np.exp(r.params["biobank"]))
    C.check("priming_or", or_ap, label="aP vs wP OR")
    C.check("priming_lci", lci, label="OR lower limit")
    C.check("priming_uci", uci, label="OR upper limit")
    C.check("priming_p", float(r.pvalues["aP"]), label="aP vs wP p")
    C.check("priming_cohort_or", or_co, label="source cohort OR")
    C.check("priming_cohort_p", float(r.pvalues["biobank"]), label="source cohort p")

    tab = np.array([[sub[sub.aP == 1]["pt_pos"].sum(), (sub.aP == 1).sum() - sub[sub.aP == 1]["pt_pos"].sum()],
                    [sub[sub.aP == 0]["pt_pos"].sum(), (sub.aP == 0).sum() - sub[sub.aP == 0]["pt_pos"].sum()]])
    or_crude, p_fisher = stats.fisher_exact(tab)

    summ["model_N"] = [len(sub)] + [""] * (len(summ) - 1)
    summ["OR_aP_vs_wP_adj"] = [round(or_ap, 3)] + [""] * (len(summ) - 1)
    summ["LCI"] = [round(lci, 3)] + [""] * (len(summ) - 1)
    summ["UCI"] = [round(uci, 3)] + [""] * (len(summ) - 1)
    summ["P_adj"] = [round(float(r.pvalues["aP"]), 4)] + [""] * (len(summ) - 1)
    summ["OR_crude"] = [round(float(or_crude), 3)] + [""] * (len(summ) - 1)
    summ["P_fisher_crude"] = [round(float(p_fisher), 4)] + [""] * (len(summ) - 1)
    summ.to_csv(C.DIR_NUM / "priming_generation.csv", index=False)
    print("\n[By generation]")
    print(summ[["generation", "birth_years", "N", "seropositive_n", "seropositive_pct",
                "GMC", "age_median", "n_main", "n_biobank"]].to_string(index=False))
    print(f"\n  main analysis (transition excluded, N={len(sub)}, adjusted for source cohort):")
    print(f"    aP vs wP  aOR {or_ap:.3f} ({lci:.3f}-{uci:.3f}) p={r.pvalues['aP']:.4f}"
          f"   | crude OR {or_crude:.3f}, Fisher p={p_fisher:.4f}")
    print(f"    source cohort OR {or_co:.3f} p={r.pvalues['biobank']:.3f}")

    # --- Sensitivity 1: two groups including the transition cohort ---
    srows = []
    for cut in C.PRIMING_SENSITIVITY_CUTS:
        for coh, s in [("pooled", p), ("main only", p[p.cohort == "main"]),
                       ("biobank only", p[p.cohort == "biobank"])]:
            a = s[s["birth_year"] >= cut]; w = s[s["birth_year"] < cut]
            if min(len(a), len(w)) < 5:
                continue
            t = np.array([[a["pt_pos"].sum(), len(a) - a["pt_pos"].sum()],
                          [w["pt_pos"].sum(), len(w) - w["pt_pos"].sum()]])
            o, pv = stats.fisher_exact(t)
            srows.append({"cutoff_birth_year": cut, "set": coh,
                          "aP_n": len(a), "aP_pct": round(100 * a["pt_pos"].mean(), 1),
                          "wP_n": len(w), "wP_pct": round(100 * w["pt_pos"].mean(), 1),
                          "OR": round(float(o), 3), "P_fisher": round(float(pv), 4)})
    sens = pd.DataFrame(srows)
    sens.to_csv(C.DIR_NUM / "priming_sensitivity.csv", index=False)
    print("\n[Sensitivity 1: two groups including the transition cohort]")
    print(sens.to_string(index=False))

    # --- Sensitivity 2: lower bound of the transition band ---
    brows = []
    for lo in C.PRIMING_TRANSITION_LOWER:
        g = p["birth_year"].map(lambda y: "aP" if y >= C.AP_BIRTH_YEAR_FROM
                                else ("transition" if y >= lo else "wP"))
        s2 = p[g.isin(["aP", "wP"])].copy()
        s2["aP"] = (g[g.isin(["aP", "wP"])] == "aP").astype(int).values
        s2["biobank"] = (s2["cohort"] == "biobank").astype(int)
        rr = sm.Logit(s2["pt_pos"], sm.add_constant(s2[["aP", "biobank"]].astype(float))).fit(disp=0)
        cc = rr.conf_int()
        brows.append({"transition_band": f"{lo}-{C.AP_BIRTH_YEAR_FROM - 1}",
                      "transition_n": int((g == "transition").sum()),
                      "wP_n": int((g == "wP").sum()), "model_N": len(s2),
                      "aOR": round(float(np.exp(rr.params["aP"])), 3),
                      "LCI": round(float(np.exp(cc.loc["aP", 0])), 3),
                      "UCI": round(float(np.exp(cc.loc["aP", 1])), 3),
                      "P": round(float(rr.pvalues["aP"]), 4)})
    band = pd.DataFrame(brows)
    band.to_csv(C.DIR_NUM / "priming_transition_band.csv", index=False)
    print("\n[Sensitivity 2: lower bound of the transition band]")
    print(band.to_string(index=False))

    # --- Figure ---
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.4),
                             gridspec_kw={"width_ratios": [1.15, 1]})
    ax = axes[0]
    x = np.arange(len(ORDER))
    vals = [summ.loc[summ.generation == g, "seropositive_pct"].iloc[0] for g in ORDER]
    ns = [summ.loc[summ.generation == g, "N"].iloc[0] for g in ORDER]
    colors = [C.PALETTE["wp"], C.PALETTE["transition"], C.PALETTE["ap"]]
    ax.bar(x, vals, 0.58, color=colors)
    for xi, v, n in zip(x, vals, ns):
        ax.text(xi, v + 1.4, f"{v:.1f}%\n(n={n})", ha="center", fontsize=9)
    ax.set_xticks(x)
    ax.set_xticklabels([LABEL[g] for g in ORDER], fontsize=9)
    ax.set_xlabel("Priming pertussis vaccine generation (by birth year)")
    ax.set_ylabel(f"Anti-PT IgG $\\geq${C.CUT_POS} {C.ASSAY_UNIT} (%)")
    ax.set_ylim(0, 88)
    ax.set_title(f"A. Seropositivity by priming generation\n"
                 f"aP vs wP: aOR {or_ap:.2f} ({lci:.2f}-{uci:.2f}), p={r.pvalues['aP']:.3f} "
                 f"(N={len(sub)}, adjusted for cohort)", fontsize=9, loc="left")
    ax.spines[["top", "right"]].set_visible(False)

    # B: same direction in each source cohort (supports pooling)
    ax = axes[1]
    w = 0.36
    for i, (coh, lab, col) in enumerate([("main", f"Clinic PWoH (n={int((p.cohort=='main').sum())})", C.PALETTE["pwoh_2025"]),
                                         ("biobank", f"Serum bank, healthy volunteers (n={int((p.cohort=='biobank').sum())})", C.PALETTE["pwoh_2025"])]):
        s = p[p["cohort"] == coh]
        v = [100 * s[s["generation"] == g]["pt_pos"].mean() for g in ORDER]
        ax.bar(np.arange(len(ORDER)) + (i - 0.5) * w, v, w, label=lab, color=col, alpha=(1.0 if coh == "main" else 0.45), hatch=(None if coh == "main" else "//"), edgecolor=col)
    ax.set_xticks(np.arange(len(ORDER)))
    ax.set_xticklabels([g if g != "transition" else "Transition" for g in ORDER], fontsize=9)
    ax.set_xlabel("Priming generation")
    ax.set_ylabel(f"Anti-PT IgG $\\geq${C.CUT_POS} {C.ASSAY_UNIT} (%)")
    ax.set_ylim(0, 88)
    ax.set_title(f"B. By source cohort\n(cohort effect OR {or_co:.2f}, p={r.pvalues['biobank']:.2f} "
                 f"— supports pooling)", fontsize=9, loc="left")
    ax.legend(frameon=False, fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)

    fig.tight_layout()
    stem = C.save_fig(fig, "priming")
    plt.close(fig)
    print(f"\n  -> numbers/priming_generation.csv, numbers/priming_sensitivity.csv, "
          f"numbers/priming_transition_band.csv, figures/{stem}.pdf")

    return {
        "priming_pooled_n": len(p), "priming_model_n": len(sub),
        "priming_ap_pct": float(summ.loc[summ.generation == "aP", "seropositive_pct"].iloc[0]),
        "priming_wp_pct": float(summ.loc[summ.generation == "wP", "seropositive_pct"].iloc[0]),
        "priming_or": round(or_ap, 3), "priming_ci": [round(lci, 3), round(uci, 3)],
        "priming_p": round(float(r.pvalues["aP"]), 4),
        "priming_cohort_p": round(float(r.pvalues["biobank"]), 3),
    }


if __name__ == "__main__":
    main()
    sys.exit(0)
