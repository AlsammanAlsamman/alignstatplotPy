"""Command-line interface: ``alignstatplot run|gui|web``."""
from __future__ import annotations

import sys

import click

from .pipeline import run_pipeline


@click.group()
@click.version_option(package_name="alignstatplot-py")
def main() -> None:
    """alignstatplot: sequence alignment statistics and fancy plots."""


@main.command()
@click.option("--fasta", "fasta_path", required=True, type=click.Path(exists=True), help="Input FASTA file.")
@click.option("--anno", "anno_path", default=None, type=click.Path(exists=True), help="Optional gene annotation file.")
@click.option("--outdir", default="output", type=click.Path(), help="Directory to write results/plots to.")
@click.option("--align-method", default="auto", help="ClustalW / ClustalOmega / Muscle / auto (informational).")
@click.option("--max-miss-per", default=0.2, type=float, help="Max fraction of gaps allowed per alignment column.")
@click.option("--n-clusters", default=4, type=int, help="Number of SNP clusters.")
@click.option("--min-cluster-length", default=3, type=int, help="Minimum SNP cluster length to draw.")
@click.option("--variant-window", default=50, type=int, help="Window size for the variant density plot.")
@click.option("--quiet", is_flag=True, help="Suppress progress logging.")
def run(fasta_path, anno_path, outdir, align_method, max_miss_per, n_clusters, min_cluster_length, variant_window, quiet):
    """Run the full analysis pipeline and save plots/tables to OUTDIR."""
    run_pipeline(
        fasta_path=fasta_path,
        anno_path=anno_path,
        outdir=outdir,
        align_method=align_method,
        max_miss_per=max_miss_per,
        n_clusters=n_clusters,
        min_cluster_length=min_cluster_length,
        variant_window=variant_window,
        verbose=not quiet,
    )
    click.echo(f"Done. Results written to {outdir}")


@main.command()
def gui() -> None:
    """Launch the desktop (Qt) GUI."""
    try:
        from .gui_qt.app import main as gui_main
    except ImportError:
        click.echo("The desktop GUI needs PySide6. Install with: pip install 'alignstatplot-py[gui]'", err=True)
        sys.exit(1)
    gui_main()


@main.command()
@click.option("--host", default="127.0.0.1")
@click.option("--port", default=8765, type=int)
@click.option("--no-browser", is_flag=True, help="Don't automatically open a browser tab.")
def web(host, port, no_browser) -> None:
    """Launch the local web GUI (FastAPI backend + browser frontend)."""
    try:
        from .webapp.server import serve
    except ImportError:
        click.echo(
            "The web GUI needs fastapi/uvicorn. Install with: pip install 'alignstatplot-py[web]'", err=True
        )
        sys.exit(1)
    serve(host=host, port=port, open_browser=not no_browser)


if __name__ == "__main__":
    main()
