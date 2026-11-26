"""
11_flowchart.py - Figure 1: participant flow and analysis sets

Every n is drawn from config.EXPECTED (frozen values). The raw row counts and the
number of second specimens removed come from the log of 00_load_and_dedup
(831 rows -> 2 second specimens of the same person removed -> 829).

Columns, left to right: PWH 2021 (purple) | PWoH (blue, two streams) |
PWH 2025 (orange) | general-population reference (grey). PWoH, the main comparison
group, sits in the centre and the external reference at the right. One box, one n.
Exclusions follow CONSORT: the main arrow runs straight down and the exclusion box
hangs to its right.

Drawn at print size (C.FIG_W wide; coordinates in inches). After drawing,
check_fit() stops the run if any text runs outside its box.

Output: figures/<FIGURES["flow"]>.pdf/.png
"""
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

import config as C

E = C.EXPECTED
RAW_2025 = E["n_2025_raw"]
DUP_2025 = E["n_2025_dup_removed"]
RAW_2021 = E["n_2021_raw"]
BB_DISPENSED = C.BIOBANK["n_dispensed"]
BB_ANALYZED = C.BIOBANK["n_analyzed"]
JIHS_ALL = E["jihs_total_n"]
JIHS_ADULT = E["jihs_pdf_n"]

INK = C.PALETTE.get("ink", "#222222")
_BOXES = []   # (patch, text) pairs, checked for overflow after drawing


def box(ax, x, y, w, h, text, fc="white", ec=INK, fs=None, bold=False, double=False):
    fs = fs or C.FS["box"]
    p = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=0.03",
                       fc=fc, ec=ec, lw=0.8)
    ax.add_patch(p)
    frame = p
    if double:
        # double outline marks the two paired-sera boxes (same 463 people)
        d = 0.035
        frame = FancyBboxPatch((x + d, y + d), w - 2 * d, h - 2 * d,
                               boxstyle="round,pad=0,rounding_size=0.02",
                               fc="none", ec=ec, lw=0.8)
        ax.add_patch(frame)
    t = ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs,
                fontweight="bold" if bold else "normal", linespacing=1.25)
    _BOXES.append((frame, t))


def arrow(ax, x0, y0, x1, y1):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                arrowprops=dict(arrowstyle="-|>", lw=0.8, color=INK, shrinkA=0, shrinkB=0,
                                mutation_scale=7))


def exclusion(ax, x_arrow, y_mid, w, h, text):
    """CONSORT: a horizontal line from the main arrow to the exclusion box on its right."""
    x_box = x_arrow + 0.12
    ax.plot([x_arrow, x_box], [y_mid, y_mid], color=INK, lw=0.8)
    box(ax, x_box, y_mid - h / 2, w, h, text, fs=C.FS["box_small"], ec="#777777")


def check_fit(fig, pad_pt=2.0):
    """Stop if any text runs outside its box (the figure is drawn at print size)."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    pad = pad_pt * fig.dpi / 72
    bad = []
    for p, t in _BOXES:
        bb, tb = p.get_window_extent(r), t.get_window_extent(r)
        if tb.x0 < bb.x0 + pad or tb.x1 > bb.x1 - pad or tb.y0 < bb.y0 + pad or tb.y1 > bb.y1 - pad:
            bad.append(t.get_text().replace("\n", " / "))
    if bad:
        raise AssertionError("text overflows its box:\n  " + "\n  ".join(bad))


def main():
    print("=" * 72)
    print("11_flowchart - Figure 1: participant flow")
    print("=" * 72)
    _BOXES.clear()
    pwh = E["n_pwh_2025"]; pwoh = E["n_pwoh_2025"]; pwoh_a = pwoh + BB_ANALYZED
    pal = C.PALETTE
    c_pwh, c_pwoh = pal["pwh_2025"], pal["pwoh_2025"]
    c_2021, c_ref = pal.get("pwh_2021", "#6A3D9A"), pal.get("jihs_2023", "#8C8C8C")
    f_pwh, f_pwoh, f_2021, f_ref = "#FBE9DD", "#E3EEF7", "#EFE6F5", "#F2F2F2"

    # Coordinates are inches at print size (C.FIG_W wide).
    # Columns, left to right: 2021 (purple) | PWoH serum bank | PWoH outpatients |
    # PWH 2025 (orange) | JIHS (grey). The exclusion boxes sit in their own row, right of
    # each main arrow, where the neighbouring column has no box.
    FW, FH = C.FIG_W, 4.3
    W, W_REF = 1.18, 1.30
    c21, cbb, cop, cpwh, cref = 0.65, 2.05, 3.35, 4.70, 6.05
    c_hosp = (cop + cpwh) / 2
    WE = {"2021": 1.00, "bb": 1.00, "hosp": 1.08, "ref": 0.95}
    Y_SRC, H_SRC = 3.42, 0.82
    Y_EXC, H_EXC = 3.00, 0.46
    Y_FORK = 2.62
    Y_GRP, H_GRP = 1.88, 0.56
    Y_SUB, H_SUB = 1.06, 0.42
    Y_SETB = 0.24

    fig = plt.figure(figsize=(FW, FH))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, FW); ax.set_ylim(0, FH); ax.axis("off")

    # ---- sources
    box(ax, c21 - W / 2, Y_SRC, W, H_SRC,
        f"PWH, previous study\nstored sera\n3 Jul 2020 – 28 Oct 2021\nn = {RAW_2021}",
        fc=f_2021, ec=c_2021)
    box(ax, cbb - W / 2, Y_SRC, W, H_SRC,
        f"IMSUT serum bank\nhealthy volunteers\n17–18 Nov 2025\nn = {BB_DISPENSED}")
    box(ax, c_hosp - W / 2, Y_SRC, W, H_SRC,
        f"IMSUT Hospital\noutpatients, residual sera\n28 Jul 2025 – 9 Jan 2026\nn = {RAW_2025}")
    box(ax, cref - W_REF / 2, Y_SRC, W_REF, H_SRC,
        f"National serosurvey\n2023 (JIHS), all ages\nn = {JIHS_ALL}", fc=f_ref, ec=c_ref)

    # ---- main arrows + exclusions
    # Neutral wording for both exclusions: the 2021 sera have no PT (or FHA) result recorded,
    # and the serum-bank sample is a laboratory assay error (Methods); the reasons differ.
    for cx, k, excl in [(c21, "2021", f"No anti-PT IgG result\nn = {RAW_2021 - E['n_pwh_2021']}"),
                        (cbb, "bb", f"No anti-PT IgG result\nn = {BB_DISPENSED - BB_ANALYZED}"),
                        (cref, "ref", f"Aged <20 years\nn = {JIHS_ALL - JIHS_ADULT}")]:
        arrow(ax, cx, Y_SRC, cx, Y_GRP + H_GRP)
        exclusion(ax, cx, Y_EXC, WE[k], H_EXC, excl)
    # hospital stream forks into PWoH outpatients and PWH 2025
    ax.plot([c_hosp, c_hosp], [Y_SRC, Y_FORK], color=INK, lw=0.8)
    exclusion(ax, c_hosp, Y_EXC, WE["hosp"], H_EXC,
              f"Second specimen of the\nsame person, n = {DUP_2025}")
    ax.plot([cop, cpwh], [Y_FORK, Y_FORK], color=INK, lw=0.8)
    arrow(ax, cop, Y_FORK, cop, Y_GRP + H_GRP)
    arrow(ax, cpwh, Y_FORK, cpwh, Y_GRP + H_GRP)

    # ---- analysis groups
    fb = C.FS["box"] + 0.5
    box(ax, c21 - W / 2, Y_GRP, W, H_GRP, f"PWH 2021 with\nanti-PT IgG\nn = {E['n_pwh_2021']}",
        fc=f_2021, ec=c_2021, fs=fb, bold=True)
    box(ax, cbb - W / 2, Y_GRP, W, H_GRP, f"PWoH serum-bank\nvolunteers\nn = {BB_ANALYZED}",
        fc=f_pwoh, ec=c_pwoh, fs=fb, bold=True)
    box(ax, cop - W / 2, Y_GRP, W, H_GRP, f"PWoH\noutpatients\nn = {pwoh}",
        fc=f_pwoh, ec=c_pwoh, fs=fb, bold=True)
    box(ax, cpwh - W / 2, Y_GRP, W, H_GRP, f"PWH 2025\nn = {pwh}",
        fc=f_pwh, ec=c_pwh, fs=fb, bold=True)
    box(ax, cref - W_REF / 2, Y_GRP, W_REF, H_GRP,
        f"General population\n≥20 years\nn = {JIHS_ADULT}", fc=f_ref, ec=c_ref, fs=fb, bold=True)

    # ---- PWoH sets span the two blue streams
    xa, wa = cbb - W / 2, (cop + W / 2) - (cbb - W / 2)
    box(ax, xa, Y_SUB, wa, H_SUB, f"PWoH set A\noutpatients + serum bank, n = {pwoh_a}",
        fc=f_pwoh, ec=c_pwoh)
    arrow(ax, cbb, Y_GRP, cbb, Y_SUB + H_SUB)
    arrow(ax, cop, Y_GRP, cop, Y_SUB + H_SUB)
    box(ax, xa, Y_SETB, wa, H_SUB, f"PWoH set B\noutpatients with blood counts, n = {pwoh}",
        fc=f_pwoh, ec=c_pwoh)
    arrow(ax, xa + wa / 2, Y_SUB, xa + wa / 2, Y_SETB + H_SUB)
    ax.text(xa + wa / 2 + 0.06, (Y_SUB + Y_SETB + H_SUB) / 2, "outpatients only",
            fontsize=C.FS["box_small"], va="center", ha="left", color="#555555")

    # ---- paired sera: the same 463 people appear under both PWH columns (2021 purple,
    # 2025 orange), at the same height. CD4/CD8 and HIV-history subsets are in Table 3.
    npair = E["mcnemar_pairs"]
    box(ax, c21 - W / 2, Y_SUB, W, H_SUB, f"Paired 2021\nn = {npair}",
        fc=f_2021, ec=c_2021, double=True)
    arrow(ax, c21, Y_GRP, c21, Y_SUB + H_SUB)
    box(ax, cpwh - W / 2, Y_SUB, W, H_SUB, f"Paired 2025\nn = {npair}",
        fc=f_pwh, ec=c_pwh, double=True)
    arrow(ax, cpwh, Y_GRP, cpwh, Y_SUB + H_SUB)

    check_fit(fig)
    stem = C.save_fig(fig, "flow")
    plt.close(fig)
    print(f"  -> figures/{stem}.pdf/.png  ({FW} x {FH} in, all text inside its box)")
    return {"flow_numbers": {"raw_2025": RAW_2025, "dup_removed": DUP_2025, "pwh": pwh,
                             "pwoh_a": pwoh_a, "pwoh_b": pwoh, "bb": BB_ANALYZED,
                             "jihs_adult": JIHS_ADULT}}


if __name__ == "__main__":
    main()
    sys.exit(0)
