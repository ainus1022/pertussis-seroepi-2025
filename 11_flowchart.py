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


def box(ax, x, y, w, h, text, fc="white", ec=INK, fs=8.5, bold=False):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.12",
                                fc=fc, ec=ec, lw=1.0))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs,
            fontweight="bold" if bold else "normal", linespacing=1.3)


def arrow(ax, x0, y0, x1, y1):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                arrowprops=dict(arrowstyle="-|>", lw=1.0, color=INK, shrinkA=0, shrinkB=0))


def exclusion(ax, x_arrow, y_mid, x_box, w, h, text):
    ax.plot([x_arrow, x_box], [y_mid, y_mid], color=INK, lw=1.0)
    box(ax, x_box, y_mid - h / 2, w, h, text, fs=7.3, ec="#777777")


def main():
    print("=" * 72)
    print("11_flowchart — Figure 1: participant flow")
    print("=" * 72)
    pwh = E["n_pwh_2025"]; pwoh = E["n_pwoh_2025"]; pwoh_a = pwoh + BB_ANALYZED
    pal = C.PALETTE
    c_pwh, c_pwoh = pal["pwh_2025"], pal["pwoh_2025"]
    c_2021, c_ref = pal.get("pwh_2021", "#6A3D9A"), pal.get("jihs_2023", "#8C8C8C")
    f_pwh, f_pwoh, f_2021, f_ref = "#FBE9DD", "#E3EEF7", "#EFE6F5", "#F2F2F2"

    # columns left -> right: 2021 (purple) | PWoH (blue, two streams) | PWH 2025 (orange) | JIHS (grey)
    W = 16.0
    X_21, X_BB, X_OP, X_PWH, X_REF = 1.0, 24.0, 43.0, 66.0, 89.0
    W_REF = 19.0
    c21, cbb, cop, cpwh, cref = (X_21 + W / 2, X_BB + W / 2, X_OP + W / 2,
                                 X_PWH + W / 2, X_REF + W_REF / 2)
    c_hosp = (cop + cpwh) / 2          # hospital source sits above the PWoH/PWH fork
    X_HOSP = c_hosp - W / 2

    fig, ax = plt.subplots(figsize=(14.4, 6.6))
    ax.set_xlim(0, 123); ax.set_ylim(26, 100); ax.axis("off")

    Y_SRC, H_SRC = 88, 11
    Y_GRP, H_GRP = 60, 11
    Y_EXC = 80.5
    WE, HE = 10.2, 7.8

    # ---- sources
    box(ax, X_21, Y_SRC, W, H_SRC,
        f"PWH, previous study\n3 Jul 2020 – 28 Oct 2021\nstored sera\nn = {RAW_2021}",
        fs=8, fc=f_2021, ec=c_2021)
    box(ax, X_BB, Y_SRC, W, H_SRC,
        f"IMSUT serum bank, 2025\nhealthy volunteers\n17–18 Nov 2025\nn = {BB_DISPENSED}", fs=8)
    box(ax, X_HOSP, Y_SRC, W, H_SRC,
        f"IMSUT Hospital\noutpatients, 2025\n28 Jul 2025 – 9 Jan 2026\nresidual sera, n = {RAW_2025}", fs=8)
    box(ax, X_REF, Y_SRC, W_REF, H_SRC,
        f"National serosurvey\n2023 (JIHS)\nall ages\nn = {JIHS_ALL}", fs=8, fc=f_ref, ec=c_ref)

    # ---- main arrows + exclusions
    y_fork = 74.0
    for cx, excl in [(c21, f"No anti-PT IgG\nresult\nn = {RAW_2021 - E['n_pwh_2021']}"),
                     (cbb, f"Not assayed\nn = {BB_DISPENSED - BB_ANALYZED}"),
                     (cref, f"Aged <20 years\nn = {JIHS_ALL - JIHS_ADULT}")]:
        arrow(ax, cx, Y_SRC, cx, Y_GRP + H_GRP)
        exclusion(ax, cx, Y_EXC, cx + W / 2 + 1.5, WE, HE, excl)
    # hospital stream forks into PWoH outpatients and PWH 2025
    ax.plot([c_hosp, c_hosp], [Y_SRC, y_fork], color=INK, lw=1.0)
    exclusion(ax, c_hosp, Y_EXC, c_hosp + W / 2 + 1.5, WE, HE,
              f"Second specimen\nof the same person\nn = {DUP_2025}")
    ax.plot([cop, cpwh], [y_fork, y_fork], color=INK, lw=1.0)
    arrow(ax, cop, y_fork, cop, Y_GRP + H_GRP)
    arrow(ax, cpwh, y_fork, cpwh, Y_GRP + H_GRP)

    # ---- analysis groups
    box(ax, X_21, Y_GRP, W, H_GRP, f"PWH 2021 with\nanti-PT IgG\nn = {E['n_pwh_2021']}",
        fc=f_2021, ec=c_2021, fs=9, bold=True)
    box(ax, X_BB, Y_GRP, W, H_GRP, f"PWoH serum-bank\nvolunteers\nn = {BB_ANALYZED}",
        fc=f_pwoh, ec=c_pwoh, fs=9, bold=True)
    box(ax, X_OP, Y_GRP, W, H_GRP, f"PWoH\noutpatients\nn = {pwoh}", fc=f_pwoh, ec=c_pwoh,
        fs=9, bold=True)
    box(ax, X_PWH, Y_GRP, W, H_GRP, f"PWH 2025\nn = {pwh}", fc=f_pwh, ec=c_pwh, fs=9, bold=True)
    box(ax, X_REF, Y_GRP, W_REF, H_GRP,
        f"General-population\ncomparator, ≥20 years\nn = {JIHS_ADULT}", fc=f_ref, ec=c_ref,
        fs=9, bold=True)

    # ---- PWoH sets span the two blue streams
    xa, wa = X_BB, (X_OP + W) - X_BB
    box(ax, xa, 46, wa, 9, f"PWoH set A\noutpatients + serum bank, n = {pwoh_a}",
        fc=f_pwoh, ec=c_pwoh, fs=8.5)
    arrow(ax, cbb, Y_GRP, cbb, 55)
    arrow(ax, cop, Y_GRP, cop, 55)
    box(ax, xa, 32, wa, 9, f"PWoH set B\noutpatients with blood counts, n = {pwoh}",
        fc=f_pwoh, ec=c_pwoh, fs=8.5)
    arrow(ax, xa + wa / 2, 46, xa + wa / 2, 41)
    ax.text(xa + wa / 2 + 1.2, 43.5, "outpatients only", fontsize=7.5, va="center",
            ha="left", color="#555555")

    # ---- PWH subsets: one box, one n
    box(ax, X_PWH, 46, W, 9, f"CD4/CD8 ratio, 2025\nn = {E['cd48_2025_n']}", fc=f_pwh, ec=c_pwh, fs=8.5)
    box(ax, X_PWH, 32, W, 9, f"HIV disease history\nn = {E['hivhx_linked_n']}", fc=f_pwh, ec=c_pwh, fs=8.5)
    arrow(ax, cpwh, Y_GRP, cpwh, 55)
    arrow(ax, cpwh, 46, cpwh, 41)

    # ---- paired sera
    box(ax, X_21, 46, W, 9, f"Paired sera\n2021 and 2025, n = {E['mcnemar_pairs']}",
        fc=f_2021, ec=c_2021, fs=8.5)
    arrow(ax, c21, Y_GRP, c21, 55)

    fig.tight_layout()
    stem = C.save_fig(fig, "flow")
    plt.close(fig)
    print(f"  -> figures/{stem}.pdf/.png")
    return {"flow_numbers": {"raw_2025": RAW_2025, "dup_removed": DUP_2025, "pwh": pwh,
                             "pwoh_a": pwoh_a, "pwoh_b": pwoh, "bb": BB_ANALYZED,
                             "jihs_adult": JIHS_ADULT}}


if __name__ == "__main__":
    main()
    sys.exit(0)
