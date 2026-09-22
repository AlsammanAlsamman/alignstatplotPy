"""Single entry point that runs the full analysis pipeline.

CLI, the Qt desktop GUI, and the web GUI all call :func:`run_pipeline`, so
results are guaranteed identical across every front end -- only R vs.
Python can differ numerically (see the README for why).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

from . import align, cluster, distance, diversity, io_utils, nucleotide, plots, stats, tree


@dataclass
class PipelineResult:
    seq_info: dict
    aligned: dict
    table: "object"
    filtered_table: "object"
    stats_df: "object"
    nuc_freq: "object"
    diversity_df: "object"
    position_stats: "object"
    region_stats_df: "object | None"
    dist_df: "object"
    tree: "object"
    cluster_result: "object"
    impact_df: "object"
    figures: dict = field(default_factory=dict)


def run_pipeline(
    fasta_path: str | Path,
    anno_path: str | Path | None = None,
    outdir: str | Path = "output",
    align_method: str = "auto",
    max_miss_per: float = 0.2,
    n_clusters: int = 4,
    min_cluster_length: int = 3,
    variant_window: int = 50,
    verbose: bool = True,
    save_figures: bool = True,
) -> PipelineResult:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    def log(msg: str) -> None:
        if verbose:
            print(f"[alignstatplot] {msg}")

    log("Reading sequences...")
    seqs = io_utils.read_fasta(fasta_path)
    anno = io_utils.read_annotation(anno_path) if anno_path else None

    log("Aligning sequences...")
    aligned = align.align_sequences(seqs, method=align_method)
    io_utils.write_fasta(aligned, outdir / "alignment.fasta")

    log("Computing per-sequence statistics...")
    stats_df = stats.per_sequence_stats(aligned)

    log("Building nucleotide table...")
    table = nucleotide.alignment_to_table(aligned)
    filtered_table = nucleotide.nuc_table_filter(table, max_miss_per=max_miss_per, remove_mono=True)
    nuc_freq = nucleotide.nuc_frequency(filtered_table)

    log("Computing diversity statistics...")
    diversity_df = diversity.diversity_stats(filtered_table)
    position_stats = diversity.position_conservation(filtered_table)
    region_stats_df = diversity.region_stats(position_stats, anno) if anno is not None else None

    log("Computing pairwise distances and tree...")
    dist_df = distance.distance_matrix(filtered_table)
    nj_tree = tree.build_tree(dist_df) if len(dist_df) > 2 else None

    log("Running SNP clustering...")
    biallelic = cluster.get_biallelic_by_freq(nuc_freq)
    ref_genotype = cluster.get_ref_genotype(nuc_freq, biallelic)
    binary_table = cluster.seq_table_to_binary(filtered_table[biallelic], ref_genotype, remove_non_ref=True)
    cluster_result = cluster.snp_cluster(binary_table, n_clusters=n_clusters) if binary_table.shape[1] > 1 else None
    impact_df = cluster.snp_cluster_impact(cluster_result) if cluster_result is not None else None

    result = PipelineResult(
        seq_info=seqs,
        aligned=aligned,
        table=table,
        filtered_table=filtered_table,
        stats_df=stats_df,
        nuc_freq=nuc_freq,
        diversity_df=diversity_df,
        position_stats=position_stats,
        region_stats_df=region_stats_df,
        dist_df=dist_df,
        tree=nj_tree,
        cluster_result=cluster_result,
        impact_df=impact_df,
    )

    log("Rendering plots...")
    figures = {
        "base_composition": plots.plot_base_composition(nuc_freq),
        "identity_distribution": plots.plot_identity_distribution(dist_df),
        "seq_stats_summary": plots.plot_seq_stats_summary(stats_df),
        "distance_heatmap": plots.distance_heatmap(dist_df),
        "conservation_track": plots.plot_conservation_track(position_stats),
        "variant_density": plots.plot_variant_density(filtered_table, window_size=variant_window),
        "sequence_logo": plots.plot_sequence_logo(filtered_table),
        "align_circle": plots.plot_align_circle(aligned),
        "summary_dashboard": plots.plot_summary_dashboard(aligned, filtered_table, dist_df, stats_df, nuc_freq),
    }
    if region_stats_df is not None and not region_stats_df.empty:
        figures["region_stats"] = plots.plot_region_stats(region_stats_df)
    if cluster_result is not None:
        figures["pca"] = plots.plot_pca(cluster_result, n_clusters=n_clusters)
        figures["snp_cluster_dendrogram"] = plots.plot_snp_cluster_dendrogram(cluster_result)
        figures["snp_cluster_pca_map"] = plots.plot_snp_cluster_pca_map(cluster_result)
    if impact_df is not None:
        elbow = cluster.find_impact_elbow(impact_df)
        figures["snp_cluster_impact"] = plots.plot_snp_cluster_impact(impact_df, elbow=elbow)

    result.figures = figures

    if save_figures:
        plot_dir = outdir / "plots"
        plot_dir.mkdir(exist_ok=True)
        for name, fig in figures.items():
            fig.savefig(plot_dir / f"{name}.png", dpi=150, bbox_inches="tight")
        stats_df.to_csv(outdir / "sequence_stats.csv", index=False)
        diversity_df.to_csv(outdir / "diversity_stats.csv", index=False)
        dist_df.to_csv(outdir / "distance_matrix.csv")
        if impact_df is not None:
            impact_df.to_csv(outdir / "snp_impact.csv", index=False)
        log(f"Saved {len(figures)} plots and result tables to {outdir}")

    return result
