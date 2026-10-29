"""
11_flowchart.py - Figure 1: participant flow and analysis sets.

All counts are taken from config.EXPECTED (frozen values). The raw row counts of the
source files (n_2025_raw, n_2025_dup_removed, n_2021_raw, jihs_total_n) are figure inputs
that match the logs of 00_load_and_dedup.py and 00b_extract_jihs.py.
Layout: four columns (2025 outpatients / serum bank / 2021 study / JIHS). Exclusions are
CONSORT-style: the main arrow runs straight down and the exclusion box hangs to its right.

Output: figures/<FIGURES["flow"]>.pdf/.png
"""
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

import config as C

E = C.EXPECTED
RAW_2025 = E["n_2025_raw"]          # 2025 records before de-duplication
DUP_2025 = E["n_2025_dup_removed"]  # second samples from the same person
RAW_2021 = E["n_2021_raw"]          # 2021 records, including those without anti-PT IgG
BB_DISPENSED = C.BIOBANK["n_dispensed"]
BB_ANALYZED = C.BIOBANK["n_analyzed"]
JIHS_ALL = E["jihs_total_n"]
JIHS_ADULT = E["jihs_pdf_n"]

INK = C.PALETTE.get("ink", "#222222")


def box(ax, x, y, w, h, text, fc="white", ec=INK, fs=8.5, bold=False):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.12",
                                fc=fc, ec=ec, lw=1.0))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs,
            fontweight="bold" if bold else "normal", linespacing=1.3)


def arrow(ax, x0, y0, x1, y1):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                arrowprops=dict(arrowstyle="-|>", lw=1.0, color=INK, shrinkA=0, shrinkB=0))


def exclusion(ax, x_arrow, y_mid, x_box, w, h, text):
    """CONSORT-style exclusion: a horizontal line from the main arrow to a box on its right."""
    ax.plot([x_arrow, x_box], [y_mid, y_mid], color=INK, lw=1.0)
    box(ax, x_box, y_mid - h / 2, w, h, text, fs=7.8, ec="#777777")


def main():
    print("=" * 72)
    print("11_flowchart — Figure 1: participant flow")
    print("=" * 72)
    pwh = E["n_pwh_2025"]; pwoh = E["n_pwoh_2025"]; pwoh_a = pwoh + BB_ANALYZED
    pal = C.PALETTE
    c_pwh, c_pwoh = pal["pwh_2025"], pal["pwoh_2025"]
    c_2021, c_ref = pal.get("pwh_2021", "#6A3D9A"), pal.get("jihs_2023", "#8C8C8C")
    f_pwh, f_pwoh, f_2021, f_ref = "#FBE9DD", "#E3EEF7", "#EFE6F5", "#F2F2F2"

    # 4 columns. In each column the source box (width WS) sits at the left,
    # the main arrow runs straight down from its centre, and the exclusion box hangs to the right.
    X = [1.0, 25.5, 50.0, 74.5]
    WS = 13.5                       # source / analysis box width
    XE = [x + WS + 2.0 for x in X]  # exclusion box x
    WE = 9.0                        # exclusion box width
    cx = [x + WS / 2 for x in X]    # arrow x

    fig, ax = plt.subplots(figsize=(12.0, 6.4))
    ax.set_xlim(0, 100); ax.set_ylim(26, 100); ax.axis("off")

    Y_SRC, H_SRC = 88, 11          # source boxes
    Y_GRP, H_GRP = 60, 11          # analysis-group boxes
    Y_EXC = 80.5                   # exclusion box centre (on the main arrow)

    # ---- sources
    box(ax, X[0], Y_SRC, WS, H_SRC, f"IMSUT Hospital\noutpatients, 2025\n28 Jul 2025 – 9 Jan 2026\nresidual sera, n = {RAW_2025}", fs=8)
    box(ax, X[1], Y_SRC, WS, H_SRC, f"IMSUT serum bank, 2025\nhealthy volunteers\n17–18 Nov 2025\nn = {BB_DISPENSED}", fs=8)
    box(ax, X[2], Y_SRC, WS, H_SRC, f"PWH, previous study\n3 Jul 2020 – 28 Oct 2021\nstored sera\nn = {RAW_2021}", fs=8, fc=f_2021, ec=c_2021)
    box(ax, X[3], Y_SRC, WS, H_SRC, f"National serosurvey\n2023 (JIHS)\nall ages\nn = {JIHS_ALL}", fs=8, fc=f_ref, ec=c_ref)

    # ---- main arrows straight down, exclusions hanging to the right
    excl = [f"Second specimen\nof the same person\nn = {DUP_2025}",
            f"Not assayed\nn = {BB_DISPENSED - BB_ANALYZED}",
            f"No anti-PT IgG\nresult\nn = {RAW_2021 - E['n_pwh_2021']}",
            f"Aged <20 years\nn = {JIHS_ALL - JIHS_ADULT}"]
    y_fork = 72.5
    # column 1 forks into PWH / PWoH outpatients below the exclusion (CONSORT-style split)
    ax.plot([cx[0], cx[0]], [Y_SRC, y_fork], color=INK, lw=1.0)
    exclusion(ax, cx[0], Y_EXC, XE[0], WE, 7.5, excl[0])
    for i in (1, 2, 3):
        arrow(ax, cx[i], Y_SRC, cx[i], Y_GRP + H_GRP)
        exclusion(ax, cx[i], Y_EXC, XE[i], WE, 7.5, excl[i])

    # ---- analysis groups
    w_half = 10.5
    x_pwh = X[0]; x_out = X[0] + w_half + 1.5
    c_pwhbox = x_pwh + w_half / 2; c_outbox = x_out + w_half / 2
    box(ax, x_pwh, Y_GRP, w_half, H_GRP, f"PWH 2025\nn = {pwh}", fc=f_pwh, ec=c_pwh, fs=9, bold=True)
    box(ax, x_out, Y_GRP, w_half, H_GRP, f"PWoH\noutpatients\nn = {pwoh}", fc=f_pwoh, ec=c_pwoh, fs=9, bold=True)
    ax.plot([c_pwhbox, c_outbox], [y_fork, y_fork], color=INK, lw=1.0)  # fork bar
    arrow(ax, c_pwhbox, y_fork, c_pwhbox, Y_GRP + H_GRP)
    arrow(ax, c_outbox, y_fork, c_outbox, Y_GRP + H_GRP)
    box(ax, X[1], Y_GRP, WS, H_GRP, f"PWoH serum-bank\nvolunteers\nn = {BB_ANALYZED}", fc=f_pwoh, ec=c_pwoh, fs=9, bold=True)
    box(ax, X[2], Y_GRP, WS, H_GRP, f"PWH 2021 with\nanti-PT IgG\nn = {E['n_pwh_2021']}", fc=f_2021, ec=c_2021, fs=9, bold=True)
    box(ax, X[3], Y_GRP, WS, H_GRP, f"General-population\ncomparator, ≥20 years\nn = {JIHS_ADULT}", fc=f_ref, ec=c_ref, fs=9, bold=True)

    # ---- PWoH sets (span outpatients + serum bank)
    xa = x_out; wa = X[1] + WS - x_out
    box(ax, xa, 44, wa, 9, f"PWoH set A (outpatients + serum bank), n = {pwoh_a}\nTable 1 · Figures 2, 3, 7", fc=f_pwoh, ec=c_pwoh, fs=8.5)
    arrow(ax, c_outbox, Y_GRP, c_outbox, 53)
    arrow(ax, cx[1], Y_GRP, cx[1], 53)
    box(ax, xa, 30, wa, 9, f"PWoH set B (outpatients with blood counts), n = {pwoh}\nTable 2 · Figure 4 (model with PWH, n = 733)", fc=f_pwoh, ec=c_pwoh, fs=8.5)
    arrow(ax, xa + wa / 2, 44, xa + wa / 2, 39)
    ax.text(xa + wa / 2 + 1.0, 41.5, "outpatients only", fontsize=7.5, va="center", ha="left", color="#555555")

    # ---- PWH subsets
    box(ax, x_pwh, 30, w_half, 23,
        f"PWH subsets\n\nCD4/CD8 ratio\n2025, n = {E['cd48_2025_n']}\n\nHIV history\nn = {E['hivhx_linked_n']}\n(Table 3)", fc=f_pwh, ec=c_pwh, fs=8)
    arrow(ax, c_pwhbox, Y_GRP, c_pwhbox, 53)

    # ---- 2021 links
    box(ax, X[2], 44, WS, 9, f"Paired sera\n2021 and 2025\nn = {E['mcnemar_pairs']} (Figure 5)", fc=f_2021, ec=c_2021, fs=8.5)
    arrow(ax, cx[2], Y_GRP, cx[2], 53)
    box(ax, X[2], 30, WS, 9, f"Baseline CD4/CD8 ratio\nlinked to 2025 titre\nn = {E['cd48_pairs']} (Figure 6)", fc=f_2021, ec=c_2021, fs=8.5)
    arrow(ax, cx[2], 44, cx[2], 39)

    fig.tight_layout()
    stem = C.save_fig(fig, "flow")
    plt.close(fig)
    print(f"  -> figures/{stem}.pdf/.png")
    return {"flow_numbers": {"raw_2025": RAW_2025, "dup_removed": DUP_2025, "pwh": pwh, "pwoh_a": pwoh_a,
                             "pwoh_b": pwoh, "bb": BB_ANALYZED, "jihs_adult": JIHS_ADULT}}


if __name__ == "__main__":
    main()
    sys.exit(0)
