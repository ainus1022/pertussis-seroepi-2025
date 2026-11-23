"""
13_titre_figures.py - Figure 8: anti-PT IgG titres

Shows the titres themselves alongside the seropositivity analyses. The reading:
exposure occurred, but titres did not rise fully (mostly subclinical infection),
and some people with HIV did not mount a rise.
Panels:
  A: PWoH set A vs PWH, individual titres (log scale) with GMC and 95% CI
  B: PWoH by priming-vaccine generation (wP / transition / aP)
  C: PWH paired sera, 2021 -> 2025
Tests on titres: Mann-Whitney U (A, B: wP vs aP) and Wilcoxon signed-rank (C).
Every GMC and P value drawn is checked against config.EXPECTED.
Output: figures/Fig8_titre.pdf / .png
"""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import config as C
import figlib as F
import lib
from matplotlib.ticker import FuncFormatter, NullFormatter

P = C.PALETTE
INK = P["ink"]
CEIL = 160.0


def gmc(v):
    v = np.asarray(v, float)
    v = v[~np.isnan(v)]
    lg = np.log(v)
    m = lg.mean(); se = lg.std(ddof=1) / np.sqrt(len(lg))
    return np.exp(m), np.exp(m - 1.96 * se), np.exp(m + 1.96 * se), len(lg)


def jitter_strip(ax, x, vals, color, width=0.17):
    rng = np.random.default_rng(20261113)
    v = np.asarray(vals, float); v = v[~np.isnan(v)]
    ax.scatter(x + rng.uniform(-width, width, len(v)), v, s=7, alpha=0.30,
               color=color, edgecolors="none", zorder=2)


def gmc_marker(ax, x, vals, color):
    g, lo, hi, n = gmc(vals)
    ax.errorbar([x], [g], yerr=[[g - lo], [hi - g]], fmt="o", ms=7, mfc="white",
                mec=color, ecolor=color, elinewidth=1.6, capsize=4, zorder=5, mew=1.8)
    # value to the right of the marker on a white chip, clear of the dots
    ax.text(x + 0.21, g, f"{g:.1f}", ha="left", va="center", fontsize=9.5,
            color=color, fontweight="bold", zorder=6,
            bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="none", alpha=0.9))
    return g, n


def refline(ax):
    for y, lab in ((10, "10"), (100, "100")):
        ax.axhline(y, color=P["ref_line"], lw=0.9, ls="--", zorder=1)
        ax.text(ax.get_xlim()[1], y, f"  {lab}", va="center", ha="left",
                fontsize=8, color=P["ref_line"])


def sig_bar(ax, x0, x1, y, label):
    ax.plot([x0, x0, x1, x1], [y, y * 1.08, y * 1.08, y], color=INK, lw=1.0)
    ax.text((x0 + x1) / 2, y * 1.13, label, ha="center", va="bottom", fontsize=9, color=INK)


def main():
    print("=" * 72); print("13_titre_figures — Figure 8: titres"); print("=" * 72)
    pwoh = lib.load_pooled_pwoh("13_titre_figures")   # via the PWoH-set guard
    assert len(pwoh) == C.EXPECTED["n_pooled_pwoh_2025"], len(pwoh)
    master = pd.read_csv(C.MASTER)
    pwh = master[master["hiv"] == 1].copy()
    pwoh["pt"] = pwoh["pt"].astype(float); pwh["pt"] = pwh["pt"].astype(float)

    key = json.loads(C.KEY_NUMBERS.read_text())

    fig, axes = plt.subplots(1, 3, figsize=(11.4, 4.3))

    # ---- A: PWoH vs PWH
    ax = axes[0]
    jitter_strip(ax, 0, pwoh["pt"], P["pwoh_2025"])
    jitter_strip(ax, 1, pwh["pt"], P["pwh_2025"])
    g0, n0 = gmc_marker(ax, 0, pwoh["pt"], P["pwoh_2025"])
    g1, n1 = gmc_marker(ax, 1, pwh["pt"], P["pwh_2025"])
    ax.set_xticks([0, 1]); ax.set_xticklabels([f"PWoH\n(n = {n0})", f"PWH\n(n = {n1})"])
    ax.set_ylabel("Anti-PT IgG (EU/mL)")
    ax.set_xlim(-0.6, 1.6)
    ax.text(0.02, 0.97, "A", transform=ax.transAxes, fontsize=13, fontweight="bold", va="top")

    # ---- B: PWoH by priming generation
    ax = axes[1]
    gens = [("wP", "wp", "wP\n(<=1974)"), ("transition", "transition", "Transition\n(1975-81)"),
            ("aP", "ap", "aP\n(>=1982)")]
    labs, gen_gmc = [], {}
    for i, (gen, ckey, lab) in enumerate(gens):
        sub = pwoh[pwoh["generation"] == gen]["pt"]
        jitter_strip(ax, i, sub, P[ckey])
        g, n = gmc_marker(ax, i, sub, P[ckey])
        gen_gmc[gen] = g
        labs.append(f"{lab}\nn = {n}")
    ax.set_xticks(range(3)); ax.set_xticklabels(labs)
    ax.set_xlim(-0.6, 2.6)
    ax.set_ylabel("Anti-PT IgG (EU/mL)")
    ax.text(0.02, 0.97, "B", transform=ax.transAxes, fontsize=13, fontweight="bold", va="top")

    # ---- C: PWH paired 2021 -> 2025
    ax = axes[2]
    pr = pwh[pwh["paired"] == 1][["pt_2021", "pt"]].astype(float).dropna()
    rng = np.random.default_rng(7)
    for _, r in pr.iterrows():
        ax.plot([0, 1], [r["pt_2021"], r["pt"]], color="#B9C2CC", lw=0.4, alpha=0.5, zorder=1)
    jitter_strip(ax, 0, pr["pt_2021"], P["pwh_2021"], width=0.10)
    jitter_strip(ax, 1, pr["pt"], P["pwh_2025"], width=0.10)
    g21, n_pair = gmc_marker(ax, 0, pr["pt_2021"], P["pwh_2021"])
    g25, _ = gmc_marker(ax, 1, pr["pt"], P["pwh_2025"])
    C.check("mcnemar_pairs", n_pair, label="Fig8C paired n")
    ax.set_xticks([0, 1]); ax.set_xticklabels([f"2021\n(n = {n_pair})", f"2025\n(n = {n_pair})"])
    ax.set_xlim(-0.5, 1.5)
    ax.set_ylabel("Anti-PT IgG (EU/mL)")
    ax.set_title("")  # no in-figure title
    ax.text(0.02, 0.97, "C", transform=ax.transAxes, fontsize=13, fontweight="bold", va="top")

    for ax in axes:
        ax.set_yscale("log")
        ax.set_ylim(0.7, 400)
        # ticks as 1 / 10 / 100 rather than powers of ten
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
        ax.yaxis.set_minor_formatter(NullFormatter())
        ax.spines[["top", "right"]].set_visible(False)
        refline(ax)

    # significance bars — titre comparisons, computed here (Mann-Whitney U on titres)
    from scipy.stats import mannwhitneyu, wilcoxon

    def pfmt(p):
        return "P < 0.001" if p < 0.001 else f"P = {p:.3f}"

    p_ab = mannwhitneyu(pwoh["pt"].dropna(), pwh["pt"].dropna()).pvalue
    p_gen = mannwhitneyu(pwoh.loc[pwoh["generation"] == "wP", "pt"].dropna(),
                         pwoh.loc[pwoh["generation"] == "aP", "pt"].dropna()).pvalue
    p_pair = wilcoxon(pr["pt_2021"], pr["pt"]).pvalue
    sig_bar(axes[0], 0, 1, 230, pfmt(p_ab))
    sig_bar(axes[1], 0, 2, 230, pfmt(p_gen))
    sig_bar(axes[2], 0, 1, 230, pfmt(p_pair))
    print(f"  titre tests: PWoH vs PWH {p_ab:.2e} / wP vs aP {p_gen:.4f} / 2021 vs 2025 paired {p_pair:.4f}")

    fig.tight_layout()
    stem = C.save_fig(fig, "titre", bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("  ->", stem)

    # ---- cross-check against frozen GMCs
    g_pwoh = gmc(pwoh["pt"])[0]; g_pwh = gmc(pwh["pt"])[0]
    print(f"  GMC PWoH set A = {g_pwoh:.2f} / PWH = {g_pwh:.2f}")
    # check every number drawn against the frozen values
    C.check_pct("pwoh_a_gmc", g_pwoh, label="Fig8A PWoH GMC")
    C.check_pct("titre_pwh_gmc", g_pwh, label="Fig8A PWH GMC")
    for gen, k in (("wP", "wp"), ("transition", "trans"), ("aP", "ap")):
        C.check_pct(f"titre_gen_{k}_gmc", gen_gmc[gen], label=f"Fig8B {gen} GMC")
    C.check_pct("titre_pair_2021_gmc", g21, label="Fig8C 2021 GMC")
    C.check_pct("titre_pair_2025_gmc", g25, label="Fig8C 2025 GMC")
    C.check("titre_mwu_p_pwoh_pwh", p_ab, label="Fig8A Mann-Whitney p")
    C.check("titre_mwu_p_wp_ap", p_gen, label="Fig8B Mann-Whitney p")
    C.check("titre_wilcoxon_p_pair", p_pair, label="Fig8C Wilcoxon p")
    print("  OK: GMCs and P values match config.EXPECTED")
    return {"titre": {"gmc_pwoh_a": round(g_pwoh, 2), "gmc_pwh": round(g_pwh, 2),
                      "gmc_gen": {k: round(v, 2) for k, v in gen_gmc.items()},
                      "gmc_pair_2021": round(g21, 2), "gmc_pair_2025": round(g25, 2), "n_pair": n_pair,
                      "p_mwu_pwoh_pwh": p_ab, "p_mwu_wp_ap": p_gen, "p_wilcoxon_pair": p_pair}}


if __name__ == "__main__":
    main()
