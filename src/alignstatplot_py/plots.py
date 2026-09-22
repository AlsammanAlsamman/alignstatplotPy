"""Matplotlib plotting functions mirroring the R package's visualisations.

Every function returns a ``matplotlib.figure.Figure`` (never draws/saves as
a side effect only), so callers -- CLI, Qt GUI, web GUI -- can all save,
embed, or further customise the same figure objects.
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, linkage

from . import circos
from .theme import ACCENT, PALETTE, apply_theme, get_seq_colors, new_figure

GAP = "-"


def plot_base_composition(nuc_freq: pd.DataFrame) -> plt.Figure:
    fig, ax = new_figure((7, 4.5))
    means = nuc_freq.mean(axis=1)
    bases = [b for b in ["A", "C", "T", "G", "N"] if b in means.index]
    colors = [PALETTE[b] for b in bases]
    ax.bar(bases, means[bases], color=colors, edgecolor="white")
    ax.set_ylabel("Mean frequency across positions")
    ax.set_title("Overall base composition")
    return fig


def plot_identity_distribution(dist_df: pd.DataFrame) -> plt.Figure:
    fig, ax = new_figure((7, 4.5))
    iu = np.triu_indices_from(dist_df.values, k=1)
    identity = 1 - dist_df.values[iu]
    ax.hist(identity, bins=20, color=ACCENT, edgecolor="white", alpha=0.85)
    ax.set_xlabel("Pairwise sequence identity")
    ax.set_ylabel("Number of pairs")
    ax.set_title("Distribution of pairwise sequence identity")
    return fig


def plot_seq_stats_summary(stats_df: pd.DataFrame) -> plt.Figure:
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    for ax in axes:
        apply_theme(ax)
    axes[0].bar(stats_df["Sequence.Name"], stats_df["Sequence.Length"], color=ACCENT)
    axes[0].set_title("Sequence length")
    axes[0].tick_params(axis="x", rotation=90)
    axes[1].bar(stats_df["Sequence.Name"], stats_df["GC.Percentage"], color=PALETTE["G"])
    axes[1].set_title("GC %")
    axes[1].tick_params(axis="x", rotation=90)
    fig.tight_layout()
    return fig


def distance_heatmap(dist_df: pd.DataFrame) -> plt.Figure:
    Z = linkage(dist_df.values[np.triu_indices_from(dist_df.values, k=1)], method="average") \
        if len(dist_df) > 2 else None
    order = list(dist_df.index)
    if Z is not None:
        dendro = dendrogram(Z, labels=order, no_plot=True)
        order = dendro["ivl"]
    ordered = dist_df.loc[order, order]

    fig, ax = plt.subplots(figsize=(6.5, 6))
    im = ax.imshow(ordered.values, cmap="viridis_r")
    ax.set_xticks(range(len(order)))
    ax.set_yticks(range(len(order)))
    fontsize = max(5, min(10, 200 / max(1, len(order))))
    ax.set_xticklabels(order, rotation=90, fontsize=fontsize)
    ax.set_yticklabels(order, fontsize=fontsize)
    fig.colorbar(im, ax=ax, label="Distance")
    ax.set_title("Pairwise distance heatmap")
    fig.tight_layout()
    return fig


def plot_conservation_track(position_stats: pd.DataFrame) -> plt.Figure:
    fig, ax = new_figure((10, 3.5))
    ax.plot(position_stats["Position"], position_stats["Conservation"], color=ACCENT, linewidth=1)
    ax.fill_between(position_stats["Position"], position_stats["Conservation"], color=ACCENT, alpha=0.2)
    ax.set_xlabel("Alignment position")
    ax.set_ylabel("Conservation")
    ax.set_ylim(0, 1.02)
    ax.set_title("Per-position conservation")
    return fig


def plot_variant_density(table: pd.DataFrame, window_size: int = 50) -> plt.Figure:
    n = len(table)
    variant = []
    for col in table.columns:
        vals = table[col]
        non_gap = vals[vals != GAP]
        variant.append(1 if non_gap.nunique() > 1 else 0)
    variant = np.array(variant)
    density = np.array(
        [variant[max(0, i - window_size // 2): i + window_size // 2].mean() for i in range(len(variant))]
    )
    fig, ax = new_figure((10, 3.5))
    ax.plot(range(1, len(density) + 1), density, color="#C62828")
    ax.set_xlabel("Alignment position")
    ax.set_ylabel(f"Variant density (window={window_size})")
    ax.set_title("Variant density along the alignment")
    return fig


def plot_region_stats(region_stats_df: pd.DataFrame) -> plt.Figure:
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    for ax in axes:
        apply_theme(ax)
    labels = region_stats_df["GeneName"].astype(str) + " " + region_stats_df["Region"].astype(str)
    axes[0].barh(labels, region_stats_df["MeanConservation"], color=ACCENT)
    axes[0].set_title("Mean conservation by region")
    axes[1].barh(labels, region_stats_df["VariantRate"], color="#C62828")
    axes[1].set_title("Variant rate by region")
    fig.tight_layout()
    return fig


def plot_pca(cluster_result, n_clusters: int = 4) -> plt.Figure:
    coords = cluster_result.pca_coords
    fig, ax = new_figure((6.5, 6))
    if coords.shape[1] >= 2:
        x, y = coords.iloc[:, 0], coords.iloc[:, 1]
    else:
        x, y = coords.iloc[:, 0], np.zeros(len(coords))
    ax.scatter(x, y, s=80, color=ACCENT, edgecolor="white", zorder=3)
    for name, xi, yi in zip(coords.index, x, y):
        ax.annotate(name, (xi, yi), fontsize=8, xytext=(4, 4), textcoords="offset points")
    var = cluster_result.explained_variance_ratio
    ax.set_xlabel(f"Dim1 ({var[0] * 100:.1f}%)" if len(var) > 0 else "Dim1")
    ax.set_ylabel(f"Dim2 ({var[1] * 100:.1f}%)" if len(var) > 1 else "Dim2")
    ax.set_title("PCA of sequences")
    return fig


def plot_snp_cluster_dendrogram(cluster_result) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(9, 5))
    apply_theme(ax)
    Z = linkage(cluster_result.loadings.values, method="ward")
    dendrogram(Z, labels=list(cluster_result.loadings.index), ax=ax, color_threshold=0)
    ax.tick_params(axis="x", labelsize=6, rotation=90)
    ax.set_title("SNP clustering (1D tree)")
    fig.tight_layout()
    return fig


def plot_snp_cluster_pca_map(cluster_result) -> plt.Figure:
    loadings = cluster_result.loadings
    clusters = cluster_result.snp_clusters
    fig, ax = new_figure((6.5, 6))
    palette = plt.cm.tab10(np.linspace(0, 1, max(clusters.nunique(), 1)))
    for i, cluster_id in enumerate(sorted(clusters.unique())):
        members = loadings.loc[clusters[clusters == cluster_id].index]
        x = members.iloc[:, 0]
        y = members.iloc[:, 1] if members.shape[1] > 1 else np.zeros(len(members))
        ax.scatter(x, y, s=40, color=palette[i], label=f"Cluster {cluster_id}", alpha=0.85)
    ax.legend(fontsize=8, frameon=False)
    ax.set_xlabel("Dim1")
    ax.set_ylabel("Dim2")
    ax.set_title("SNP cluster PCA map")
    return fig


def plot_snp_cluster_impact(impact_df: pd.DataFrame, elbow: int | None = None) -> plt.Figure:
    sorted_df = impact_df.sort_values("Impact", ascending=False).reset_index(drop=True)
    cumulative = sorted_df["Impact"].cumsum() / sorted_df["Impact"].sum()

    fig, ax1 = new_figure((8, 4.5))
    ax1.bar(range(1, len(sorted_df) + 1), sorted_df["Impact"], color=ACCENT, alpha=0.7)
    ax1.set_xlabel("SNP rank")
    ax1.set_ylabel("Impact", color=ACCENT)
    ax2 = ax1.twinx()
    ax2.plot(range(1, len(sorted_df) + 1), cumulative, color="#C62828")
    ax2.set_ylabel("Cumulative impact share", color="#C62828")
    if elbow:
        ax1.axvline(elbow, color="black", linestyle="--", linewidth=1)
    ax1.set_title("SNP cluster impact")
    fig.tight_layout()
    return fig


def plot_sequence_logo(table: pd.DataFrame, max_positions: int = 60) -> plt.Figure:
    """A simple, dependency-free sequence logo: letters stacked and scaled
    by per-position frequency."""
    from .nucleotide import nuc_frequency

    freq = nuc_frequency(table)
    positions = list(freq.columns)[:max_positions]
    fig, ax = plt.subplots(figsize=(max(6, len(positions) * 0.25), 4))
    apply_theme(ax)
    for xi, pos in enumerate(positions):
        col = freq[pos].sort_values()
        y0 = 0.0
        for base, height in col.items():
            if height <= 0:
                continue
            ax.text(
                xi,
                y0 + height / 2,
                base,
                fontsize=max(6, min(22, 300 * height)),
                ha="center",
                va="center",
                color=PALETTE.get(base, "#333333"),
                fontweight="bold",
            )
            y0 += height
    ax.set_xlim(-1, len(positions))
    ax.set_ylim(0, 1)
    ax.set_xticks(range(len(positions)))
    ax.set_xticklabels(positions, rotation=90, fontsize=6)
    ax.set_ylabel("Frequency")
    ax.set_title("Sequence logo")
    fig.tight_layout()
    return fig


def plot_align_circle(
    aligned: dict[str, str],
    colors: list[str] | None = None,
    max_gene_sectors: int = 15,
) -> plt.Figure:
    """Circular alignment overview. Python analogue of the R package's
    ``plotAlignCircle``: dispatches to :func:`draw_cons_with_genes` (one
    pie sector per sequence plus a consensus sector, linked by ribbons)
    for up to ``max_gene_sectors`` sequences -- matching R's own 15-sequence
    threshold, above which the per-sector ribbon layout gets too crowded to
    read -- and to :func:`draw_cons_with_no_genes` (stacked concentric
    coverage rings around a single consensus circle) beyond that.
    """
    if len(aligned) <= max_gene_sectors:
        return draw_cons_with_genes(aligned, colors=colors)
    return draw_cons_with_no_genes(aligned, colors=colors)


def draw_cons_with_genes(
    aligned: dict[str, str],
    cons_zoom_factor: float = 3.0,
    colors: list[str] | None = None,
    seq_fontsize: float = 8,
    cons_fontsize: float = 10,
    tick_fontsize: float = 5,
    cons_tick_fontsize: float = 6,
    link_alpha: float = 0.4,
    max_fragments_per_seq: int = 40,
) -> plt.Figure:
    """One pie sector per sequence (sized by its ungapped length) plus a
    zoomed-up ``Consensus`` sector, connected by translucent ribbon links
    at every ungapped fragment -- and a matching stack of opaque per-sequence
    coverage rings inside the consensus sector. Python analogue of the R
    package's ``drawConsWithGenes``, built from the same fragment/link
    algorithm (``alignmentNoGaps`` / ``alignmentNoGapsLinks``) rather than
    ``circlize``."""
    names = list(aligned.keys())
    n = len(names)
    cons_length = len(next(iter(aligned.values())))
    seq_lengths = {name: max(len(seq.replace(GAP, "")), 1) for name, seq in aligned.items()}

    sizes = {**seq_lengths, "Consensus": cons_length * cons_zoom_factor}
    order = names + ["Consensus"]
    sectors = circos.layout_sectors(sizes, order, gap_degree=5.0, start_degree=90.0)
    # `sizes["Consensus"]` is zoomed only to *allocate more angular width* to
    # the consensus sector -- its actual coordinate range for tick/ribbon
    # placement is still 0..cons_length, not 0..cons_length*zoom.
    zoomed = sectors["Consensus"]
    cons_sector = circos.Sector("Consensus", zoomed.start_deg, zoomed.end_deg, 0.0, cons_length)
    sectors["Consensus"] = cons_sector
    palette = colors or get_seq_colors(n)

    fig, ax = plt.subplots(figsize=(9, 9))
    ax.set_xlim(-1.35, 1.35)
    ax.set_ylim(-1.35, 1.35)
    ax.set_aspect("equal")
    ax.axis("off")

    ring_band = 0.5 / max(n, 1)
    link_inner_radius = 0.39
    ring_outer_max = 1 - 0.11

    for i, name in enumerate(names):
        sector = sectors[name]
        color = palette[i % len(palette)]

        fragments = circos.non_gap_fragments(aligned[name])
        local_fragments = circos.local_fragment_positions(fragments)
        frag_pairs = circos.merge_fragment_pairs(list(zip(fragments, local_fragments)), max_gap=2)
        frag_pairs = circos.limit_fragment_pairs(frag_pairs, max_n=max_fragments_per_seq)

        ring_outer = ring_outer_max - ring_band * i
        ring_inner = ring_outer_max - ring_band * (i + 1)

        for (a_start, a_end), (l_start, l_end) in frag_pairs:
            cons_start = cons_length - a_end
            cons_end = cons_length - a_start + 1

            circos.draw_ribbon(
                ax, cons_sector, cons_start, cons_end, link_inner_radius, sector, l_start, l_end, 0.9,
                facecolor=color, edgecolor="none", alpha=link_alpha, zorder=1,
            )
            circos.draw_arc_band(
                ax, cons_sector, cons_start, cons_end, ring_inner, ring_outer,
                facecolor=color, edgecolor="white", linewidth=0.3, zorder=3,
            )

        circos.draw_sector_ticks(ax, sector, 0.92, n_ticks=3, fontsize=tick_fontsize)
        circos.draw_sector_label(ax, sector, 1.06, name, fontsize=seq_fontsize, fontweight="bold", color="#222222")

    circos.draw_sector_ticks(ax, cons_sector, 0.92, n_ticks=6, fontsize=cons_tick_fontsize, fmt="{:.0f}bp")
    circos.draw_sector_label(
        ax, cons_sector, 1.06, "Consensus", fontsize=cons_fontsize, fontweight="bold", color="#333333"
    )

    ax.set_title("Circular alignment overview", pad=14, fontsize=13, fontweight="bold")
    return fig


def draw_cons_with_no_genes(
    aligned: dict[str, str],
    colors: list[str] | None = None,
    bg_color: str = "#CCCCCC",
    seq_fontsize: float | None = None,
    bp_fontsize: float = 7,
    max_fragments_per_seq: int = 60,
) -> plt.Figure:
    """A single consensus circle (one small angular gap) with one stacked
    coverage ring per sequence -- each ring highlighting that sequence's
    ungapped stretches in its own colour over a light background baseline.
    Python analogue of the R package's ``drawConsWithNoGenes``, used once a
    per-sequence pie-sector layout would get too crowded to read."""
    names = list(aligned.keys())
    n = len(names)
    cons_length = len(next(iter(aligned.values())))
    palette = colors or get_seq_colors(n)
    if seq_fontsize is None:
        seq_fontsize = max(4.0, min(11.0, 90.0 / max(n, 1)))

    sectors = circos.layout_sectors({"Consensus": cons_length}, ["Consensus"], gap_degree=15.0, start_degree=95.0)
    cons_sector = sectors["Consensus"]

    fig, ax = plt.subplots(figsize=(9, 9))
    ax.set_xlim(-1.15, 1.15)
    ax.set_ylim(-1.15, 1.15)
    ax.set_aspect("equal")
    ax.axis("off")

    r_min, r_max = 0.18, 0.95
    band = (r_max - r_min) / max(n, 1)

    for i, name in enumerate(names):
        r_inner = r_min + band * i
        r_outer = r_min + band * (i + 1)
        color = palette[i % len(palette)]

        circos.draw_arc_band(
            ax, cons_sector, 0, cons_length, r_inner, r_outer, facecolor=bg_color, edgecolor="none", zorder=1
        )

        fragments = circos.merge_close_fragments(circos.non_gap_fragments(aligned[name]), max_gap=2)
        if len(fragments) > max_fragments_per_seq:
            fragments = sorted(fragments, key=lambda f: f[1] - f[0], reverse=True)[:max_fragments_per_seq]
            fragments.sort()
        for a_start, a_end in fragments:
            circos.draw_arc_band(
                ax, cons_sector, a_start, a_end, r_inner, r_outer,
                facecolor=color, edgecolor="white", linewidth=0.15, zorder=2,
            )

        theta = cons_sector.angle_of(cons_sector.data_min)
        mid_r = (r_inner + r_outer) / 2
        x, y = mid_r * np.cos(theta), mid_r * np.sin(theta)
        deg = np.rad2deg(theta)
        # Run the label tangentially (perpendicular to the radius) rather
        # than radially, so stacked rings' labels sit in their own row
        # instead of overlapping along the same radial line.
        tangent = deg - 90
        flipped = tangent < -90 or tangent > 90
        ax.text(
            x, y, name, fontsize=seq_fontsize, ha=("left" if flipped else "right"), va="center",
            rotation=tangent + 180 if flipped else tangent, rotation_mode="anchor",
        )

    circos.draw_sector_ticks(ax, cons_sector, r_max + 0.02, n_ticks=8, fontsize=bp_fontsize, fmt="{:.0f}bp")
    ax.set_title("Circular alignment overview", pad=14, fontsize=13, fontweight="bold")
    return fig


def plot_summary_dashboard(aligned, table, dist_df, stats_df, nuc_freq) -> plt.Figure:
    """A single composite figure combining several panels (R's
    ``plotSummaryDashboard``, built with ``patchwork`` there)."""
    fig = plt.figure(figsize=(14, 9))
    gs = fig.add_gridspec(2, 3)

    ax1 = fig.add_subplot(gs[0, 0])
    apply_theme(ax1)
    means = nuc_freq.mean(axis=1)
    bases = [b for b in ["A", "C", "T", "G", "N"] if b in means.index]
    ax1.bar(bases, means[bases], color=[PALETTE[b] for b in bases])
    ax1.set_title("Base composition")

    ax2 = fig.add_subplot(gs[0, 1])
    apply_theme(ax2)
    iu = np.triu_indices_from(dist_df.values, k=1)
    ax2.hist(1 - dist_df.values[iu], bins=15, color=ACCENT)
    ax2.set_title("Identity distribution")

    ax3 = fig.add_subplot(gs[0, 2])
    apply_theme(ax3)
    ax3.bar(stats_df["Sequence.Name"], stats_df["GC.Percentage"], color=PALETTE["G"])
    ax3.tick_params(axis="x", rotation=90, labelsize=6)
    ax3.set_title("GC %")

    ax4 = fig.add_subplot(gs[1, :2])
    im = ax4.imshow(dist_df.values, cmap="viridis_r")
    ax4.set_xticks(range(len(dist_df)))
    ax4.set_yticks(range(len(dist_df)))
    ax4.set_xticklabels(dist_df.index, rotation=90, fontsize=6)
    ax4.set_yticklabels(dist_df.index, fontsize=6)
    ax4.set_title("Distance heatmap")
    fig.colorbar(im, ax=ax4, fraction=0.046)

    ax5 = fig.add_subplot(gs[1, 2])
    apply_theme(ax5)
    ax5.bar(stats_df["Sequence.Name"], stats_df["Sequence.Length"], color=ACCENT)
    ax5.tick_params(axis="x", rotation=90, labelsize=6)
    ax5.set_title("Sequence length")

    fig.suptitle("alignstatplot summary dashboard", fontsize=14, fontweight="bold")
    fig.tight_layout()
    return fig
