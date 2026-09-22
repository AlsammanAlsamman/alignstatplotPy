"""Shared visual theme for all alignstatplot-py figures (mirrors R's
``theme_alignstatplot``)."""
from __future__ import annotations

import matplotlib.pyplot as plt
from matplotlib.colors import to_hex

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


def get_seq_colors(n: int) -> list[str]:
    """A distinct categorical colour per sequence, analogous to the R
    package's ``getSeqColors`` (Dark2 / Kelly / Alphabet / Polychrome by
    sequence count)."""
    if n <= 8:
        cmap = plt.get_cmap("Dark2")
        return [to_hex(cmap(i)) for i in range(max(n, 1))]
    if n <= 20:
        cmap = plt.get_cmap("tab20")
        return [to_hex(cmap(i)) for i in range(n)]
    pool: list[str] = []
    for name in ("tab20", "tab20b", "tab20c"):
        cmap = plt.get_cmap(name)
        pool += [to_hex(cmap(i)) for i in range(20)]
    if n <= len(pool):
        return pool[:n]
    reps = n // len(pool) + 1
    return (pool * reps)[:n]
