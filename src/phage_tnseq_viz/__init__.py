"""phage-tnseq-viz: Tn-Seq processing hand-off and phage genome visualisation."""

from __future__ import annotations

__version__ = "0.2.0"

from .dataset import FinalSite, load_final_dataset
from .essentiality import InsertionSite, annotate_sites_with_genes, classify_genes
from .genome import Gene, GenomeRecord, gene_identifier, load_genome
from .pipeline import PipelineConfig, ReadPair, ToolPaths, run_pipeline
from .plot import PlotOptions, render
from .transposon import Transposon, find_insertion_sites, resolve_transposon

__all__ = [
    "__version__",
    "Gene",
    "GenomeRecord",
    "gene_identifier",
    "load_genome",
    "FinalSite",
    "load_final_dataset",
    "InsertionSite",
    "annotate_sites_with_genes",
    "classify_genes",
    "PipelineConfig",
    "ReadPair",
    "ToolPaths",
    "run_pipeline",
    "PlotOptions",
    "render",
    "Transposon",
    "resolve_transposon",
    "find_insertion_sites",
]
