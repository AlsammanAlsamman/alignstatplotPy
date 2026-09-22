"""Shared visual theme for all alignstatplot-py figures (mirrors R's
``theme_alignstatplot``)."""
from __future__ import annotations

import matplotlib.pyplot as plt

PALETTE = {
    "A": "#2E7D32",
    "C": "#1565C0",
    "T": "#C62828",
    "G": "#F9A825",
    "N": "#9E9E9E",
    "-": "#E0E0E0",
}
ACCENT = "#5E35B1"
BACKGROUND = "#FAFAFA"


def apply_theme(ax) -> None:
    ax.set_facecolor(BACKGROUND)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.tick_params(colors="#333333")
    ax.grid(axis="y", color="#E0E0E0", linewidth=0.6, zorder=0)
    ax.set_axisbelow(True)


def new_figure(figsize=(8, 5)):
    fig, ax = plt.subplots(figsize=figsize)
    fig.patch.set_facecolor("white")
    apply_theme(ax)
    return fig, ax
