"""A small circos-style rendering engine (sectors, ribbons, arcs, radial
ticks) -- just enough of what the R package leans on from ``circlize`` to
reproduce its two alignment-circle layouts in Matplotlib, with no
Bioconductor/circlize equivalent needed in Python.
"""
from __future__ import annotations

import re
from dataclasses import dataclass

import numpy as np
from matplotlib.patches import PathPatch, Wedge
from matplotlib.path import Path

GAP_CHAR = "-"


def non_gap_fragments(aligned_seq: str) -> list[tuple[int, int]]:
    """1-based inclusive (start, end) runs of non-gap characters, in
    alignment (gapped) coordinates -- mirrors R's ``alignmentNoGaps``."""
    return [(m.start() + 1, m.end()) for m in re.finditer(r"[^-]+", aligned_seq)]


def local_fragment_positions(fragments: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Map each alignment-coordinate fragment to its position range in the
    *ungapped* sequence -- mirrors R's ``alignmentNoGapsLinks``. Fragments
    are contiguous non-gap runs in alignment order, so ungapped positions
    accumulate directly."""
    local = []
    cursor = 1
    for start, end in fragments:
        length = end - start + 1
        local.append((cursor, cursor + length - 1))
        cursor += length
    return local


@dataclass
class Sector:
    name: str
    start_deg: float  # inclusive
    end_deg: float  # inclusive, start_deg >= end_deg (clockwise layout)
    data_min: float
    data_max: float

    def angle_of(self, value: float) -> float:
        span = self.data_max - self.data_min
        frac = 0.0 if span == 0 else (value - self.data_min) / span
        frac = min(max(frac, 0.0), 1.0)
        deg = self.start_deg + frac * (self.end_deg - self.start_deg)
        return np.deg2rad(deg)


def layout_sectors(
    sizes: "dict[str, float]", order: list[str], gap_degree: float = 5.0, start_degree: float = 90.0
) -> "dict[str, Sector]":
    """Lay sectors clockwise around the circle, sized proportional to
    ``sizes``, separated by a fixed ``gap_degree``."""
    n = len(order)
    available = 360.0 - n * gap_degree
    total = sum(sizes[name] for name in order)
    sectors: "dict[str, Sector]" = {}
    cursor = start_degree
    for name in order:
        span = (sizes[name] / total) * available if total else available / n
        sectors[name] = Sector(name, cursor, cursor - span, 0.0, sizes[name])
        cursor -= span + gap_degree
    return sectors


def _xy(theta: np.ndarray, r: float) -> np.ndarray:
    return np.column_stack([r * np.cos(theta), r * np.sin(theta)])


def draw_arc_band(ax, sector: Sector, value_start: float, value_end: float, r_inner: float, r_outer: float, **kwargs):
    """A filled annular band over ``[value_start, value_end]`` of a sector."""
    theta1 = np.rad2deg(sector.angle_of(value_start))
    theta2 = np.rad2deg(sector.angle_of(value_end))
    lo, hi = sorted((theta1, theta2))
    wedge = Wedge((0, 0), r_outer, lo, hi, width=r_outer - r_inner, **kwargs)
    ax.add_patch(wedge)
    return wedge


def draw_ribbon(
    ax,
    sector_a: Sector,
    a_start: float,
    a_end: float,
    radius_a: float,
    sector_b: Sector,
    b_start: float,
    b_end: float,
    radius_b: float,
    n: int = 30,
    **kwargs,
):
    """A filled ribbon connecting an arc on ``sector_a`` to an arc on
    ``sector_b``, curving through the centre -- the Matplotlib analogue of
    ``circlize``'s ``circos.genomicLink``."""
    a1, a2 = sector_a.angle_of(a_start), sector_a.angle_of(a_end)
    b1, b2 = sector_b.angle_of(b_start), sector_b.angle_of(b_end)

    arc_a = _xy(np.linspace(a1, a2, n), radius_a)
    arc_b = _xy(np.linspace(b1, b2, n), radius_b)

    verts = [arc_a[0]]
    codes = [Path.MOVETO]
    for pt in arc_a[1:]:
        verts.append(pt)
        codes.append(Path.LINETO)
    verts += [(0.0, 0.0), arc_b[0]]
    codes += [Path.CURVE3, Path.CURVE3]
    for pt in arc_b[1:]:
        verts.append(pt)
        codes.append(Path.LINETO)
    verts += [(0.0, 0.0), arc_a[0]]
    codes += [Path.CURVE3, Path.CURVE3]
    verts.append(arc_a[0])
    codes.append(Path.CLOSEPOLY)

    patch = PathPatch(Path(verts, codes), **kwargs)
    ax.add_patch(patch)
    return patch


def draw_sector_ticks(ax, sector: Sector, radius: float, n_ticks: int = 5, fontsize: float = 5, fmt="{:.0f}"):
    """Small radial tick marks + labels at the outer edge of a sector,
    each rotated to stay legible ("niceFacing")."""
    values = np.linspace(sector.data_min, sector.data_max, n_ticks)
    for v in values:
        theta = sector.angle_of(v)
        x0, y0 = radius * np.cos(theta), radius * np.sin(theta)
        x1, y1 = (radius + 0.015) * np.cos(theta), (radius + 0.015) * np.sin(theta)
        ax.plot([x0, x1], [y0, y1], color="#333333", linewidth=0.7, zorder=5)
        deg = np.rad2deg(theta)
        rotation = deg if -90 <= deg <= 90 else deg + 180
        ha = "left" if -90 <= deg <= 90 else "right"
        lx, ly = (radius + 0.03) * np.cos(theta), (radius + 0.03) * np.sin(theta)
        ax.text(
            lx, ly, fmt.format(v), fontsize=fontsize, rotation=rotation, rotation_mode="anchor",
            ha=ha, va="center", zorder=5,
        )


def merge_close_fragments(fragments: list[tuple[int, int]], max_gap: int = 2) -> list[tuple[int, int]]:
    """Merge fragments separated by a short gap (<= ``max_gap``) into one,
    purely for rendering simplicity (fewer, cleaner arcs)."""
    if not fragments:
        return []
    merged = [fragments[0]]
    for start, end in fragments[1:]:
        last_start, last_end = merged[-1]
        if start - last_end - 1 <= max_gap:
            merged[-1] = (last_start, end)
        else:
            merged.append((start, end))
    return merged


def merge_fragment_pairs(
    frag_pairs: "list[tuple[tuple[int, int], tuple[int, int]]]", max_gap: int = 2
) -> "list[tuple[tuple[int, int], tuple[int, int]]]":
    """Like :func:`merge_close_fragments`, but keeps a parallel (aligned
    coordinate, local ungapped coordinate) pair in lock-step, for the
    ribbon-link layout where both coordinate systems matter."""
    if not frag_pairs:
        return []
    merged = [frag_pairs[0]]
    for (a_start, a_end), (l_start, l_end) in frag_pairs[1:]:
        (pa_start, pa_end), (pl_start, pl_end) = merged[-1]
        if a_start - pa_end - 1 <= max_gap:
            merged[-1] = ((pa_start, a_end), (pl_start, l_end))
        else:
            merged.append(((a_start, a_end), (l_start, l_end)))
    return merged


def limit_fragment_pairs(frag_pairs, max_n: int = 40):
    """Keep only the ``max_n`` longest fragments (by aligned span), to
    bound rendering cost on messy/heavily-gapped alignments."""
    if len(frag_pairs) <= max_n:
        return frag_pairs
    ranked = sorted(frag_pairs, key=lambda p: p[0][1] - p[0][0], reverse=True)[:max_n]
    return sorted(ranked, key=lambda p: p[0][0])


def draw_sector_label(ax, sector: Sector, radius: float, label: str, fontsize: float = 8, **kwargs):
    mid = (sector.start_deg + sector.end_deg) / 2
    theta = np.deg2rad(mid)
    x, y = radius * np.cos(theta), radius * np.sin(theta)
    rotation = mid if -90 <= mid <= 90 else mid + 180
    ha = "left" if -90 <= mid <= 90 else "right"
    ax.text(x, y, label, fontsize=fontsize, rotation=rotation, rotation_mode="anchor", ha=ha, va="center", **kwargs)
