"""Local web GUI: a FastAPI backend calling the same
:func:`alignstatplot_py.pipeline.run_pipeline` used by the CLI and Qt GUI,
plus a static HTML/CSS/JS frontend (see ``frontend/`` in the repo root for
its source) for uploading FASTA files via drag-and-drop and browsing
results in the browser.
"""
from __future__ import annotations

import base64
import io
import tempfile
import threading
import webbrowser
from pathlib import Path

import pandas as pd
from fastapi import FastAPI, File, Form, UploadFile
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from fastapi.staticfiles import StaticFiles

from ..pipeline import run_pipeline

STATIC_DIR = Path(__file__).parent / "static"

app = FastAPI(title="alignstatplot")
app.add_middleware(
    CORSMiddleware, allow_origins=["*"], allow_methods=["*"], allow_headers=["*"]
)


def _fig_to_base64(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=140, bbox_inches="tight")
    return base64.b64encode(buf.getvalue()).decode("ascii")


def _df_to_records(df: pd.DataFrame) -> list[dict]:
    return df.reset_index().to_dict(orient="records")


@app.post("/api/run")
async def run_analysis(
    fasta: UploadFile = File(...),
    annotation: UploadFile | None = File(None),
    max_miss_per: float = Form(0.2),
    n_clusters: int = Form(4),
    min_cluster_length: int = Form(3),
    variant_window: int = Form(50),
):
    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        fasta_path = tmp_path / (fasta.filename or "input.fasta")
        fasta_path.write_bytes(await fasta.read())

        anno_path = None
        if annotation is not None and annotation.filename:
            anno_path = tmp_path / annotation.filename
            anno_path.write_bytes(await annotation.read())

        try:
            result = run_pipeline(
                fasta_path=fasta_path,
                anno_path=anno_path,
                outdir=tmp_path / "output",
                max_miss_per=max_miss_per,
                n_clusters=n_clusters,
                min_cluster_length=min_cluster_length,
                variant_window=variant_window,
                verbose=False,
                save_figures=False,
            )
        except Exception as exc:  # surfaced to the frontend as a normal error message
            return JSONResponse(status_code=400, content={"error": str(exc)})

        plots = {name: _fig_to_base64(fig) for name, fig in result.figures.items()}
        return {
            "stats": _df_to_records(result.stats_df),
            "diversity": _df_to_records(result.diversity_df),
            "distance": _df_to_records(result.dist_df),
            "impact": _df_to_records(result.impact_df) if result.impact_df is not None else [],
            "plots": plots,
        }


if STATIC_DIR.exists():
    app.mount("/", StaticFiles(directory=STATIC_DIR, html=True), name="static")


def serve(host: str = "127.0.0.1", port: int = 8765, open_browser: bool = True) -> None:
    import uvicorn

    url = f"http://{host}:{port}"
    if open_browser:
        threading.Timer(1.0, lambda: webbrowser.open(url)).start()
    print(f"[alignstatplot] Web GUI running at {url}")
    uvicorn.run(app, host=host, port=port, log_level="warning")


if __name__ == "__main__":
    serve()
