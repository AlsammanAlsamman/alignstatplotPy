from pathlib import Path

import pandas as pd

from alignstatplot_py import io_utils, nucleotide, stats
from alignstatplot_py.pipeline import run_pipeline

DATA_DIR = Path(__file__).parent / "data"


def test_read_fasta_already_aligned():
    seqs = io_utils.read_fasta(DATA_DIR / "Example_Sequences_Aligned.fasta")
    assert len(seqs) > 0
    lengths = {len(s) for s in seqs.values()}
    assert len(lengths) == 1  # pre-aligned fixture is equal-length


def test_per_sequence_stats_columns():
    seqs = io_utils.read_fasta(DATA_DIR / "Example_Sequences_Aligned.fasta")
    df = stats.per_sequence_stats(seqs)
    for col in ["SeqNo", "Sequence.Name", "Sequence.Length", "A", "T", "C", "G", "Gap", "GC.Percentage"]:
        assert col in df.columns
    assert len(df) == len(seqs)


def test_nuc_table_filter_drops_gappy_and_mono_columns():
    seqs = io_utils.read_fasta(DATA_DIR / "Example_Sequences_Aligned.fasta")
    table = nucleotide.alignment_to_table(seqs)
    filtered = nucleotide.nuc_table_filter(table, max_miss_per=0.2, remove_mono=True)
    assert filtered.shape[1] <= table.shape[1]
    # every remaining column must be polymorphic among non-gap bases
    for col in filtered.columns:
        non_gap = filtered[col][filtered[col] != "-"]
        assert non_gap.nunique() > 1


def test_run_pipeline_end_to_end(tmp_path):
    result = run_pipeline(
        fasta_path=DATA_DIR / "Example_Small.fasta",
        anno_path=DATA_DIR / "Example_Small.anno",
        outdir=tmp_path,
        n_clusters=3,
        verbose=False,
    )
    assert isinstance(result.stats_df, pd.DataFrame)
    assert len(result.stats_df) == len(result.seq_info)
    assert not result.diversity_df.empty
    assert result.dist_df.shape[0] == result.dist_df.shape[1]
    assert "base_composition" in result.figures
    assert "align_circle" in result.figures
    assert (tmp_path / "plots").exists()
    assert (tmp_path / "sequence_stats.csv").exists()
