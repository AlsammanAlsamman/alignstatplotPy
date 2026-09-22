"""Binary genotype conversion and PCA-based SNP clustering (mirrors R's
``seqTableToBinary``, ``SNPCluster``, and ``SNPClusterImpact``)."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from sklearn.decomposition import PCA

GAP = "-"
VALID_BASES = ["A", "C", "T", "G"]


def get_biallelic_by_freq(nuc_freq: pd.DataFrame) -> list[str]:
    """Positions where exactly two of A/C/T/G are observed."""
    bases = nuc_freq.loc[[b for b in VALID_BASES if b in nuc_freq.index]]
    return [col for col in bases.columns if (bases[col] > 0).sum() == 2]


def get_ref_genotype(nuc_freq: pd.DataFrame, biallelic: list[str]) -> dict[str, str]:
    """The majority allele at each biallelic position, used as the
    reference for binary (0/1) recoding."""
    bases = nuc_freq.loc[[b for b in VALID_BASES if b in nuc_freq.index]]
    return {col: bases[col].idxmax() for col in biallelic}


def seq_table_to_binary(table: pd.DataFrame, ref, remove_non_ref: bool = True) -> pd.DataFrame:
    """Recode a nucleotide table to 0 (matches reference)/1 (alternate).

    ``ref`` is either a dict of ``{position: reference_allele}`` (typical:
    from :func:`get_ref_genotype`) or a list of sequence names to use as a
    reference genotype (one or more references; positions shared/agreeing
    across all reference sequences supply the reference allele).
    """
    if isinstance(ref, dict):
        columns = [c for c in table.columns if c in ref]
        out = pd.DataFrame(index=table.index, columns=columns, dtype=object)
        for col in columns:
            ref_allele = ref[col]
            for name in table.index:
                val = table.loc[name, col]
                if val == GAP:
                    out.loc[name, col] = GAP
                else:
                    out.loc[name, col] = 0 if val == ref_allele else 1
    else:
        ref_names = list(ref)
        ref_rows = table.loc[ref_names]
        columns = []
        ref_allele_map = {}
        for col in table.columns:
            alleles = set(ref_rows[col]) - {GAP}
            if len(alleles) == 1:
                columns.append(col)
                ref_allele_map[col] = alleles.pop()
        out = pd.DataFrame(index=table.index, columns=columns, dtype=object)
        for col in columns:
            ref_allele = ref_allele_map[col]
            for name in table.index:
                val = table.loc[name, col]
                out.loc[name, col] = GAP if val == GAP else (0 if val == ref_allele else 1)

    if remove_non_ref:
        keep = [c for c in out.columns if set(out[c]) - {GAP} <= {0, 1}]
        out = out[keep]
    return out


@dataclass
class ClusterResult:
    binary_table: pd.DataFrame
    pca_coords: pd.DataFrame
    loadings: pd.DataFrame
    explained_variance_ratio: np.ndarray
    snp_clusters: pd.Series
    model: PCA = field(repr=False)


def snp_cluster(binary_table: pd.DataFrame, n_clusters: int = 4, n_components: int = 5) -> ClusterResult:
    """PCA over SNPs (columns), then hierarchical clustering of SNPs in PCA
    (loading) space -- the Python analogue of R's MCA-based ``SNPCluster``."""
    numeric = binary_table.replace(GAP, np.nan).astype(float)
    numeric = numeric.fillna(numeric.mean())

    n_components = max(1, min(n_components, numeric.shape[0] - 1, numeric.shape[1]))
    pca = PCA(n_components=n_components)
    sample_coords = pca.fit_transform(numeric.T.values)  # SNPs as observations
    loadings = pd.DataFrame(
        sample_coords, index=numeric.columns, columns=[f"Dim{i + 1}" for i in range(n_components)]
    )

    n_clusters = max(1, min(n_clusters, len(loadings) - 1)) if len(loadings) > 1 else 1
    if n_clusters > 1:
        Z = linkage(loadings.values, method="ward")
        cluster_ids = fcluster(Z, t=n_clusters, criterion="maxclust")
    else:
        cluster_ids = np.ones(len(loadings), dtype=int)
    clusters = pd.Series(cluster_ids, index=loadings.index, name="Cluster")

    seq_pca = PCA(n_components=min(n_components, numeric.shape[0] - 1) or 1)
    seq_coords = seq_pca.fit_transform(numeric.values)
    seq_coords_df = pd.DataFrame(
        seq_coords, index=numeric.index, columns=[f"Dim{i + 1}" for i in range(seq_coords.shape[1])]
    )

    return ClusterResult(
        binary_table=binary_table,
        pca_coords=seq_coords_df,
        loadings=loadings,
        explained_variance_ratio=pca.explained_variance_ratio_,
        snp_clusters=clusters,
        model=pca,
    )


def snp_cluster_impact(result: ClusterResult) -> pd.DataFrame:
    """Score every SNP by its contribution to the retained PCA dimensions,
    weighted by how much variance each dimension explains."""
    loadings = result.loadings
    weights = result.explained_variance_ratio[: loadings.shape[1]]
    weights = weights / weights.sum()

    sq = loadings.values**2
    row_sums = sq.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1
    cos2 = (sq / row_sums) @ weights
    impact = sq @ weights

    df = pd.DataFrame(
        {"SNP": loadings.index, "Impact": impact, "Cos2": cos2, "Cluster": result.snp_clusters.values}
    )
    return df.sort_values("Impact", ascending=False).reset_index(drop=True)


def find_impact_elbow(impact_df: pd.DataFrame) -> int:
    """Index (1-based count of SNPs to keep) at the elbow of the sorted
    cumulative-impact curve, via the max-distance-to-chord heuristic."""
    values = impact_df.sort_values("Impact", ascending=False)["Impact"].to_numpy()
    n = len(values)
    if n <= 2:
        return n
    x = np.arange(n)
    y = values
    p1 = np.array([x[0], y[0]])
    p2 = np.array([x[-1], y[-1]])
    line_vec = p2 - p1
    line_len = np.linalg.norm(line_vec)
    line_unit = line_vec / line_len if line_len else line_vec
    distances = []
    for xi, yi in zip(x, y):
        p = np.array([xi, yi]) - p1
        proj_len = np.dot(p, line_unit)
        proj = proj_len * line_unit
        distances.append(np.linalg.norm(p - proj))
    return int(np.argmax(distances)) + 1


def filter_high_impact_snps(
    table: pd.DataFrame, impact_df: pd.DataFrame, top_n: int | None = None, min_impact: float | None = None
) -> pd.DataFrame:
    if top_n is not None:
        keep = impact_df.sort_values("Impact", ascending=False).head(top_n)["SNP"]
    elif min_impact is not None:
        keep = impact_df[impact_df["Impact"] >= min_impact]["SNP"]
    else:
        elbow = find_impact_elbow(impact_df)
        keep = impact_df.sort_values("Impact", ascending=False).head(elbow)["SNP"]
    keep = [c for c in table.columns if c in set(keep)]
    return table[keep].copy()
