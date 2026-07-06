"""phage-tnseq-viz: pre-visualization for phage Tn-Seq (HIDDEN-seq) experiments."""

from __future__ import annotations

__version__ = "0.1.0"

from .genome import Gene, GenomeRecord, load_genome
from .plot import PlotOptions, render
from .transposon import Transposon, find_insertion_sites, resolve_transposon

__all__ = [
    "__version__",
    "Gene",
    "GenomeRecord",
    "load_genome",
    "PlotOptions",
    "render",
    "Transposon",
    "resolve_transposon",
    "find_insertion_sites",
]
