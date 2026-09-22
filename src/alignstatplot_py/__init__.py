"""alignstatplot-py: sequence alignment statistics and fancy plots.

Python port/companion of the R package ``alignstatplot``
(https://github.com/AlsammanAlsamman/alignstatplot). See ``pipeline.run_pipeline``
for the single entry point used by the CLI and both GUIs.
"""
from .pipeline import PipelineResult, run_pipeline

__version__ = "0.1.0"
__all__ = ["run_pipeline", "PipelineResult", "__version__"]
