"""Matplotlib plotting functions mirroring the R package's visualisations.

Every function returns a ``matplotlib.figure.Figure`` (never draws/saves as
a side effect only), so callers -- CLI, Qt GUI, web GUI -- can all save,
embed, or further customise the same figure objects.
"""
from __future__ import annotations

import math

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, Rectangle, Wedge
from scipy.cluster.hierarchy import dendrogram, linkage

from .theme import ACCENT, PALETTE, apply_theme, new_figure

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


def plot_align_circle(aligned: dict[str, str], annotation: pd.DataFrame | None = None) -> plt.Figure:
    """Circular alignment plot: each sequence is a concentric ring coloured
    by base, with an optional gene-annotation ring at the centre. Python
    analogue of the R package's ``plotAlignCircle`` / ``drawConsWithGenes``."""
    names = list(aligned.keys())
    n = len(names)
    length = len(next(iter(aligned.values())))

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection="polar")
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_ylim(0, n + 2)
    ax.axis("off")

    theta = np.linspace(0, 2 * np.pi, length, endpoint=False)
    width = (2 * np.pi / length) * 1.02

    for ring, name in enumerate(reversed(names)):
        seq = aligned[name]
        radius = ring + 1
        colors = [PALETTE.get(b, "#333333") for b in seq]
        ax.bar(theta, height=0.9, width=width, bottom=radius, color=colors, linewidth=0)
        ax.text(0, radius + 0.45, name, fontsize=6, ha="right", va="center")

    if annotation is not None and not annotation.empty:
        for _, region in annotation.iterrows():
            start_theta = 2 * np.pi * (region["Start"] - 1) / length
            end_theta = 2 * np.pi * (region["End"] - 1) / length
            color = "#5E35B1" if str(region.get("Region", "")).lower() == "exon" else "#B39DDB"
            ax.bar(
                (start_theta + end_theta) / 2,
                height=0.6,
                width=max(end_theta - start_theta, width),
                bottom=0.2,
                color=color,
                linewidth=0,
            )

    ax.set_title("Circular alignment overview", pad=20)
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
