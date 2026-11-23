"""
figlib.py - shared figure helpers

  panel()         no in-figure title; only the panel letter (A/B) at the top left
  sig_bar()       significance bracket with its P value (linear or log axis)
  colour_legend() legend text coloured like its series
  p_label()       one P-value format for every figure
"""
import numpy as np
from matplotlib.ticker import FuncFormatter, NullFormatter
import config as C

INK = C.PALETTE.get("ink", "#222222")


def panel(ax, letter, fs=None):
    """Remove the in-figure title and place the panel letter only."""
    fs = fs or C.FS["panel"]
    ax.set_title("")
    ax.text(-0.02, 1.04, letter, transform=ax.transAxes, fontsize=fs,
            fontweight="bold", va="bottom", ha="right", color=INK)


def p_label(p):
    p = float(p)
    if p < 0.001:
        return "P < 0.001"
    return f"P = {p:.3f}" if p < 0.1 else f"P = {p:.2f}"


def sig_bar(ax, x1, x2, y, p=None, text=None, fs=None, color=INK):
    """Draw a bracket over x1..x2 labelled with the P value (or text)."""
    fs = fs or C.FS["annot"]
    label = text if text is not None else p_label(p)
    if ax.get_yscale() == "log":
        h = y * 0.10
        ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=0.8, color=color, clip_on=False)
        ax.text((x1 + x2) / 2, (y + h) * 1.04, label, ha="center", va="bottom",
                fontsize=fs, color=color, clip_on=False)
    else:
        lo, hi = ax.get_ylim()
        h = (hi - lo) * 0.020
        ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=0.8, color=color, clip_on=False)
        ax.text((x1 + x2) / 2, y + h * 1.25, label, ha="center", va="bottom",
                fontsize=fs, color=color, clip_on=False)


def _handle_colour(h):
    for attr in ("get_facecolor", "get_color"):
        f = getattr(h, attr, None)
        if f is None:
            continue
        try:
            c = f()
        except Exception:
            continue
        if c is None:
            continue
        a = np.asarray(c, dtype=float).reshape(-1)
        if a.size >= 3:
            return (float(a[0]), float(a[1]), float(a[2]), 1.0)   # opaque, so the text stays readable
    return None


def colour_legend(ax, **kw):
    """Draw the legend and colour each label like its series."""
    kw.setdefault("frameon", False)
    leg = ax.legend(**kw)
    for txt, h in zip(leg.get_texts(), leg.legend_handles):
        col = _handle_colour(h)
        if col is not None:
            txt.set_color(col)
    return leg


def log_ticks(ax, axis="y"):
    """Label a log axis 1 / 10 / 100 instead of powers of ten."""
    a = ax.yaxis if axis == "y" else ax.xaxis
    a.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    a.set_minor_formatter(NullFormatter())
