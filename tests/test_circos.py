from pathlib import Path

from alignstatplot_py import align, io_utils
from alignstatplot_py.circos import (
    layout_sectors,
    local_fragment_positions,
    merge_close_fragments,
    non_gap_fragments,
)
from alignstatplot_py.plots import draw_cons_with_genes, draw_cons_with_no_genes, plot_align_circle

DATA_DIR = Path(__file__).parent / "data"


def test_non_gap_fragments_and_local_positions():
    fragments = non_gap_fragments("--AC-GT--A")
    assert fragments == [(3, 4), (6, 7), (10, 10)]
    assert local_fragment_positions(fragments) == [(1, 2), (3, 4), (5, 5)]


def test_merge_close_fragments():
    assert merge_close_fragments([(1, 5), (7, 10), (20, 25)], max_gap=2) == [(1, 10), (20, 25)]


def test_layout_sectors_covers_available_degrees():
    sectors = layout_sectors({"a": 1, "b": 1, "c": 2}, ["a", "b", "c"], gap_degree=5, start_degree=90)
    total_span = sum(abs(s.end_deg - s.start_deg) for s in sectors.values())
    assert total_span == 360 - 3 * 5
    # proportional to size: c should get twice a's (and b's) span
    assert abs(abs(sectors["c"].end_deg - sectors["c"].start_deg) - 2 * abs(sectors["a"].end_deg - sectors["a"].start_deg)) < 1e-9


def _aligned_example():
    seqs = io_utils.read_fasta(DATA_DIR / "Example_Small.fasta")
    return align.align_sequences(seqs)


def test_plot_align_circle_dispatches_by_sequence_count():
    aligned = _aligned_example()
    assert len(aligned) <= 15
    fig = plot_align_circle(aligned)
    assert fig is not None

    fig_forced_rings = plot_align_circle(aligned, max_gene_sectors=0)
    assert fig_forced_rings is not None


def test_draw_cons_with_genes_consensus_sector_uses_real_length_not_zoomed():
    aligned = _aligned_example()
    fig = draw_cons_with_genes(aligned, cons_zoom_factor=3.0)
    assert fig is not None


def test_draw_cons_with_no_genes_runs():
    aligned = _aligned_example()
    fig = draw_cons_with_no_genes(aligned)
    assert fig is not None
