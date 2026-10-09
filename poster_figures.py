"""
poster_figures.py - figures for the IDWeek 2026 poster (P-1187).

These are separate renderings of the same results as the manuscript figures,
laid out for a poster: no in-figure titles, one colour per population
(config.PALETTE), and the three main bar charts (1A, 2A, 3A) share figure size,
axes rectangle, y-axis and physical bar width.

Values are read from numbers/*.csv. Only panel 3B needs participant-level data
(the private analysis file).

Output: figures/poster/Poster_Fig1A ... Poster_Fig3B (PNG at 300 dpi, and PDF)
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
import config as C

OUT = C.DIR_FIG / "poster"
OUT.mkdir(exist_ok=True)
P = C.PALETTE

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": 20, "axes.labelsize": 20, "xtick.labelsize": 19, "ytick.labelsize": 19,
    "legend.fontsize": 19, "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 1.2,
})

# The three main bar charts use a fixed axes rectangle (no tight_layout) so that
# axis width and height are identical across them.
TRIO_FIGSIZE = (9.5, 5.8)
TRIO_AXRECT = [0.105, 0.185, 0.875, 0.705]
TRIO_BARFRAC = 0.055      # width of one bar as a fraction of the axis width


def trio_axes():
    fig = plt.figure(figsize=TRIO_FIGSIZE)
    ax = fig.add_axes(TRIO_AXRECT)
    return fig, ax


def trio_grouped(ax, cats, v1, v2, c1, c2, l1, l2):
    n = len(cats)
    w = TRIO_BARFRAC * n      # x-range is n wide, so this keeps the physical width constant
    x = np.arange(n)
    b1 = ax.bar(x - w / 2, v1, w, color=c1, label=l1)
    b2 = ax.bar(x + w / 2, v2, w, color=c2, label=l2)
    for b in list(b1) + list(b2):      # value labels inside the bars
        ax.text(b.get_x() + b.get_width() / 2, b.get_height() - 1.5, f"{b.get_height():.0f}",
                ha="center", va="top", fontsize=17, color="white", fontweight="bold")
    ax.set_xlim(-0.6, n - 0.4)
    ax.set_xticks(x); ax.set_xticklabels(cats)
    ax.set_ylim(0, 90); ax.set_yticks(range(0, 81, 20))   # headroom for significance bars
    ax.set_ylabel(f"Anti-PT IgG ≥{C.CUT_POS} {C.ASSAY_UNIT} (%)")
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), frameon=False, ncol=2)
    return x, w, b1, b2


def save(fig, name, tight=False):
    if tight:
        fig.tight_layout()
    fig.savefig(OUT / f"{name}.png", dpi=300, facecolor="white")
    fig.savefig(OUT / f"{name}.pdf", facecolor="white")
    plt.close(fig)
    print("  ->", (OUT / f"{name}.png").relative_to(C.REPO))


def stars(p):
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""


# ---------- 1A: general population 2023 vs PWoH 2025, by age ----------
def fig1a():
    d = pd.read_csv(C.DIR_NUM / "analysis1_agestrat.csv")
    fig, ax = trio_axes()
    cats = [f"{g}\n({j}/{p})" for g, j, p in zip(d.age_grp, d.JIHS_N, d.PWoH_N)]
    trio_grouped(ax, cats, d.JIHS_pct, d.PWoH_pct, P["jihs_2023"], P["pwoh_2025"],
                 "General population, JIHS 2023", "PWoH 2025 (this study)")
    ax.set_xlabel("Age group, years  (JIHS n / PWoH n)")
    save(fig, "Poster_Fig1A_jihs_vs_pwoh_by_age")


# ---------- 1B: PWoH by priming-vaccine generation ----------
def fig1b():
    d = pd.read_csv(C.DIR_NUM / "priming_generation.csv").set_index("generation")
    order = ["wP", "transition", "aP"]
    labs = ["wP\n(≤1974)", "Transition\n(1975–81)", "aP\n(≥1982)"]
    cols = [P["pwoh_2025"]] * 3      # one population, one colour
    fig = plt.figure(figsize=(5.6, 5.8))
    ax = fig.add_axes([0.20, 0.185, 0.77, 0.75])
    bars = ax.bar(range(3), d.loc[order, "seropositive_pct"], 0.64, color=cols)
    for b, g in zip(bars, order):
        ax.text(b.get_x() + b.get_width() / 2, b.get_height() - 1.5,
                f"{b.get_height():.1f}%\n(n={int(d.loc[g,'N'])})", ha="center", va="top",
                fontsize=17, color="white", fontweight="bold")
    ax.set_xticks(range(3)); ax.set_xticklabels(labs)
    ax.set_xlim(-0.6, 2.6)
    ax.set_ylim(0, 90); ax.set_yticks(range(0, 81, 20))
    ax.set_ylabel(f"Anti-PT IgG ≥{C.CUT_POS} {C.ASSAY_UNIT} (%)")
    ax.set_xlabel("Priming vaccine (birth year)")
    top = max(d.loc[["wP", "aP"], "seropositive_pct"]) + 3.0
    ax.plot([0, 0, 2, 2], [top, top + 2.0, top + 2.0, top], color="black", lw=1.4)
    ax.text(1, top + 2.1, stars(float(d.loc["wP", "P_adj"])), ha="center", va="bottom",
            fontsize=22, fontweight="bold")
    save(fig, "Poster_Fig1B_priming_generation")


# ---------- 2A: PWoH vs PWH by age, with significance bars ----------
def fig2a():
    d = pd.read_csv(C.DIR_NUM / "table1_agestrat.csv")
    fig, ax = trio_axes()
    cats = [f"{g}\n({a}/{b})" for g, a, b in zip(d.age_grp, d.PWoH_N, d.PWH_N)]
    x, w, b1, b2 = trio_grouped(ax, cats, d.PWoH_pct, d.PWH_pct, P["pwoh_2025"], P["pwh_2025"],
                                f"PWoH 2025 (n={int(d.PWoH_N.sum())})",
                                f"PWH 2025 (n={int(d.PWH_N.sum())})")
    for i, p in enumerate(d.P_fisher):
        s = stars(p)
        if s:
            top = max(d.PWoH_pct[i], d.PWH_pct[i]) + 2.0
            ax.plot([x[i] - w / 2, x[i] - w / 2, x[i] + w / 2, x[i] + w / 2],
                    [top, top + 1.5, top + 1.5, top], color="black", lw=1.4)
            ax.text(x[i], top + 1.6, s, ha="center", va="bottom", fontsize=20, fontweight="bold")
    ax.set_xlabel("Age group, years  (PWoH n / PWH n)")
    save(fig, "Poster_Fig2A_pwoh_vs_pwh_by_age")


# ---------- 2B: forest plot ----------
def fig2b():
    d = pd.read_csv(C.DIR_NUM / "analysis2.csv")
    labels = ["HIV-positive\n(PWH vs PWoH)", "Age\n(per 10 years)", "Male sex", "Lymphocyte count\n(per 500 /µL)"]
    fig = plt.figure(figsize=(8.8, 5.8))
    ax = fig.add_axes([0.32, 0.20, 0.36, 0.72])
    y = np.arange(len(d))[::-1]
    for yi, (_, r) in zip(y, d.iterrows()):
        hiv = r.Covariate.startswith("HIV")
        col = P["pwh_2025"] if hiv else P["ink"]
        ax.errorbar(r.aOR, yi, xerr=[[r.aOR - r.LCI], [r.UCI - r.aOR]], fmt="o", color=col, ecolor=col,
                    capsize=5, markersize=11, lw=2.2)
        ax.annotate(f"{r.aOR:.2f} ({r.LCI:.2f}–{r.UCI:.2f})\np={r.P:.3f}", xy=(1.06, yi),
                    xycoords=("axes fraction", "data"), va="center", ha="left", fontsize=19,
                    fontweight="bold" if hiv else "normal", color=col if hiv else "black")
    ax.axvline(1.0, color=P["ref_line"], ls="--", lw=1.5)
    ax.set_xscale("log"); ax.set_xlim(0.3, 2.6); ax.set_xticks([0.5, 1.0, 2.0]); ax.set_xticklabels(["0.5", "1.0", "2.0"])
    ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.set_ylim(-0.6, len(d) - 0.4)
    ax.set_yticks(y); ax.set_yticklabels(labels)
    ax.set_xlabel(f"Adjusted OR for anti-PT IgG ≥{C.CUT_POS}\n{C.ASSAY_UNIT} (95% CI)")
    ax.spines["left"].set_visible(False); ax.tick_params(axis="y", length=0)
    save(fig, "Poster_Fig2B_forest")


# ---------- 3A: PWH paired 2021 vs 2025, by age ----------
def fig3a():
    d = pd.read_csv(C.DIR_NUM / "mcnemar_agegroup.csv")
    fig, ax = trio_axes()
    cats = [f"{g}\n(n={n})" for g, n in zip(d.age_grp, d.n)]
    trio_grouped(ax, cats, d.pos_2021_pct, d.pos_2025_pct, P["pwh_2021"], P["pwh_2025"], "PWH 2021", "PWH 2025")
    ax.set_xlabel("Age group in 2025, years")
    save(fig, "Poster_Fig3A_pwh_paired_by_age")


# ---------- 3B: 2025 titre by baseline CD4/CD8 ratio ----------
def fig3b():
    m = pd.read_csv(C.MASTER)
    p = m[(m["hiv"] == 1) & (m["cd48_linked"] == 1)].copy()
    lo = p[p["cd48_ratio_2021"] <= C.CD48_CUT]["pt"]
    hi = p[p["cd48_ratio_2021"] > C.CD48_CUT]["pt"]
    fig, ax = plt.subplots(figsize=(6.8, 5.6))
    bp = ax.boxplot([lo, hi], positions=[0, 1], widths=0.55, patch_artist=True, showfliers=False,
                    medianprops=dict(color=P["ink"], lw=2.2))
    for patch, a in zip(bp["boxes"], [0.30, 0.65]):
        patch.set_facecolor(P["pwh_2025"]); patch.set_alpha(a); patch.set_edgecolor(P["ink"])
    rng = np.random.default_rng(0)
    for i, s in enumerate([lo, hi]):
        ax.scatter(rng.normal(i, 0.06, len(s)), s, s=12, alpha=0.35, color=P["ink"], lw=0)
    ax.axhline(C.CUT_POS, color=P["ref_line"], ls="--", lw=1.5)
    ax.set_yscale("log"); ax.set_ylim(0.8, 250)
    ax.set_yticks([1, 10, 100]); ax.set_yticklabels(["1", "10", "100"]); ax.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.set_xticks([0, 1]); ax.set_xticklabels([f"≤0.5\n(n={len(lo)})", f">0.5\n(n={len(hi)})"])
    ax.set_xlabel("Baseline CD4/CD8 ratio (2021)")
    ax.set_ylabel(f"Anti-PT IgG in 2025 ({C.ASSAY_UNIT}, log)")
    # significance bar for low vs high stratum (numbers/cd4cd8.csv)
    pv = float(pd.read_csv(C.DIR_NUM / "cd4cd8.csv").loc[0, "P_mwu"])
    top = 150
    ax.plot([0, 0, 1, 1], [top, top * 1.18, top * 1.18, top], color="black", lw=1.4)
    ax.text(0.5, top * 1.20, stars(pv), ha="center", va="bottom", fontsize=24, fontweight="bold")
    save(fig, "Poster_Fig3B_cd4cd8_titre", tight=True)


if __name__ == "__main__":
    print("poster_figures ->", OUT.relative_to(C.REPO))
    fig1a(); fig1b(); fig2a(); fig2b(); fig3a(); fig3b()
    print("done")
