from pathlib import Path

from alignstatplot_py import cluster, distance, diversity, io_utils, nucleotide

DATA_DIR = Path(__file__).parent / "data"


def _filtered_table():
    seqs = io_utils.read_fasta(DATA_DIR / "Example_Sequences_Aligned.fasta")
    table = nucleotide.alignment_to_table(seqs)
    return nucleotide.nuc_table_filter(table, max_miss_per=0.2, remove_mono=True)


def test_diversity_stats_sane_ranges():
    table = _filtered_table()
    df = diversity.diversity_stats(table)
    row = df.iloc[0]
    assert row["NumSequences"] == len(table)
    assert 0 <= row["NucleotideDiversity"] <= 1
    assert row["SegregatingSites"] <= row["AlignmentLength"]


def test_distance_matrix_symmetric_zero_diagonal():
    table = _filtered_table()
    dist = distance.distance_matrix(table)
    assert (dist.values.diagonal() == 0).all()
    assert (dist.values == dist.values.T).all()


def test_snp_cluster_and_impact():
    table = _filtered_table()
    freq = nucleotide.nuc_frequency(table)
    biallelic = cluster.get_biallelic_by_freq(freq)
    assert len(biallelic) > 0
    ref = cluster.get_ref_genotype(freq, biallelic)
    binary = cluster.seq_table_to_binary(table[biallelic], ref, remove_non_ref=True)
    assert set(binary.stack().unique()) <= {0, 1, "-"}

    result = cluster.snp_cluster(binary, n_clusters=3)
    assert result.snp_clusters.nunique() <= 3
    impact = cluster.snp_cluster_impact(result)
    assert list(impact["SNP"]) == list(impact.sort_values("Impact", ascending=False)["SNP"])

    filtered = cluster.filter_high_impact_snps(table[biallelic], impact, top_n=5)
    assert filtered.shape[1] == 5
