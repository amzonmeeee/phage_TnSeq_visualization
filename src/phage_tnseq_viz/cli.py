"""Command-line interface for phage-tnseq-viz.

The original one-command pre-visualisation interface remains valid.  Two
explicit paths extend it for sequencing data:

* ``plot`` (or the legacy bare command) accepts a user-supplied final-site CSV;
* ``process`` optionally runs FastQC, fastp, TRANSIT TPP+BWA, and SeqKit before
  writing that final CSV and rendering it.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, replace
import json
from pathlib import Path
import shlex
import subprocess
import sys
from typing import Iterable, Mapping

from . import __version__
from .dataset import (
    DatasetError,
    FinalSite,
    fill_missing_final_sites,
    final_sites_from_counts,
    group_counts_for_plotting,
    load_final_dataset,
    write_final_dataset,
)
from .essentiality import (
    AnnotatedInsertionSite,
    ClassificationResult,
    InsertionSite,
    annotate_sites_with_genes,
    classify_genes,
    load_classifier,
)
from .genome import GenomeRecord, load_genome
from .pipeline import (
    PipelineConfig,
    PipelineError,
    ReadPair,
    ToolPaths,
    candidate_insertion_sites,
    run_pipeline,
)
from .plot import PAPER_SIZES, PlotOptions, render
from .transposon import PRESETS, Transposon, find_insertion_sites, resolve_transposon


def build_parser() -> argparse.ArgumentParser:
    """Build the backwards-compatible bare ``phage-tnseq-viz GENOME.gbk`` parser."""
    p = argparse.ArgumentParser(
        prog="phage-tnseq-viz",
        description=(
            "Draw a phage genome map, optionally overlaying a processed Tn-Seq "
            "final dataset. For raw FASTQ processing run `phage-tnseq-viz process --help`."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    _add_plot_arguments(p)
    return p


def build_plot_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="phage-tnseq-viz plot",
        description="Plot a GenBank reference with optional final Tn-Seq counts.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    _add_plot_arguments(p)
    return p


def build_process_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="phage-tnseq-viz process",
        description=(
            "Optionally process Illumina Tn-Seq FASTQ reads with FastQC, fastp, "
            "TRANSIT TPP+BWA and SeqKit, then write final CSV/map outputs."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("input", help="Reference phage genome in GenBank (.gbk/.gb) format.")
    p.add_argument("--reads1", required=True, help="Read 1 FASTQ(.gz).")
    p.add_argument("--reads2", default=None, help="Optional read 2 FASTQ(.gz).")
    p.add_argument("-o", "--output-dir", default="tnseq_run", help="Run-output directory.")
    p.add_argument("--sample-id", default="sample", help="Stable label used in intermediate file names.")
    p.add_argument("--interactive", action="store_true", help="Prompt for the main optional stages and parameters.")

    stages = p.add_argument_group("optional processing stages")
    stages.add_argument("--skip-fastqc", action="store_true", help="Do not run FastQC initial QC.")
    stages.add_argument("--skip-fastp", action="store_true", help="Do not run fastp trimming/filtering.")
    stages.add_argument("--skip-tpp", action="store_true", help="Do not run TPP+BWA (no final count table is produced).")
    stages.add_argument("--quality-phred", type=int, default=20, help="fastp qualified-base Phred cutoff.")
    stages.add_argument("--min-read-length", type=int, default=20, help="Discard reads shorter than this after fastp.")
    stages.add_argument("--adapter-sequence", default=None, help="Known R1 adapter sequence for fastp.")
    stages.add_argument("--adapter-sequence-r2", default=None, help="Known R2 adapter sequence for fastp.")
    stages.add_argument("--threads", type=int, default=1, help="FastQC/fastp worker threads.")

    tpp = p.add_argument_group("TPP + BWA junction mapping")
    tpp.add_argument("--tpp-bin", "--tpp", dest="tpp_bin", default="tpp.py", help="TPP executable or tpp.py path.")
    tpp.add_argument("--tpp-python", default=None, help="Interpreter to prepend when --tpp-bin is a tpp.py script.")
    tpp.add_argument("--bwa-bin", default="bwa", help="BWA executable passed to TPP.")
    tpp.add_argument("--fastqc-bin", default="fastqc", help="FastQC executable.")
    tpp.add_argument("--fastp-bin", default="fastp", help="fastp executable.")
    tpp.add_argument("--seqkit-bin", default="seqkit", help="SeqKit executable.")
    tpp.add_argument("--tpp-mode", choices=("himar1", "tn5"), default="himar1", help="TPP protocol/library family.")
    tpp.add_argument("--tpp-primer", default=None, help="Actual transposon terminal prefix to find in R1; required for non-standard IR cargo libraries.")
    tpp.add_argument("--tpp-mismatches", type=int, default=1, help="Allowed prefix/constant-region mismatches in TPP.")
    tpp.add_argument("--tpp-replicon-ids", default=None, help="Comma-separated TPP output contig IDs for multi-contig references.")
    tpp.add_argument("--min-mapped-reads", "--read-count-threshold", type=float, default=0, help="Set sites with fewer mapped reads to zero.")

    subsampling = p.add_argument_group("optional depth matching with SeqKit")
    subsampling.add_argument("--subsample-depth", type=int, default=None, help="Randomly sample this many reads before each TPP run.")
    subsampling.add_argument("--subsample-replicates", type=int, default=1, help="Number of seeded depth-matched subsamples to average.")
    subsampling.add_argument("--subsample-seed", type=int, default=11, help="First deterministic SeqKit random seed.")
    subsampling.add_argument("--observed-sites-only", action="store_true", help="Do not add zero-count candidate sites (not recommended for essentiality).")

    final = p.add_argument_group("final data and map")
    final.add_argument("--classifier", default=None, help="Trusted local .py custom classifier exposing classify_gene(gene_id, site_rows).")
    final.add_argument("--no-essentiality", action="store_true", help="Write/plot counts but skip gene essentiality calls.")
    final.add_argument("--no-read-histogram", action="store_true", help="Do not draw blue per-site read-count bars.")
    final.add_argument("--read-histogram-cap", type=float, default=None, metavar="PCT", help="Cap the read-histogram scale at this percentile (e.g. 95) of each contig's positive counts; taller bars clip to full height and the scale reads '≥ N'.")
    final.add_argument("--no-plot", action="store_true", help="Write final CSV(s) but do not render an image.")
    final.add_argument("--plot-output", default="tnseq_map.png", help="Map filename, relative to --output-dir unless absolute.")
    final.add_argument("--show-candidate-sites", action="store_true", help="Also draw potential TA insertion ticks on the processed-data map.")
    return p


def _add_plot_arguments(p: argparse.ArgumentParser) -> None:
    p.add_argument("input", help="Input phage genome in GenBank (.gbk/.gb) format.")
    p.add_argument("-o", "--output", default="phage_map.png", help="Output image path (.png or .svg).")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    g_tn = p.add_argument_group("transposon / candidate sites")
    g_tn.add_argument(
        "-t", "--transposon", default="mariner",
        help="Preset name (%s) or custom IUPAC insertion motif (e.g. TA, NTAN)."
        % ", ".join(sorted(PRESETS)),
    )
    g_tn.add_argument("--no-insertion-sites", action="store_true", help="Do not draw theoretical insertion-site ticks.")
    g_tn.add_argument("--show-theoretical-sites", action="store_true", help="With --final-dataset, also show legacy theoretical-site ticks/density.")
    g_tn.add_argument("--no-insertion-density", action="store_true", help="Do not draw the theoretical insertion-density heat track.")
    g_tn.add_argument("--density-window", type=int, default=None, help="Density window in bp; default is genome length/100.")
    g_tn.add_argument("--single-strand", action="store_true", help="Scan a custom motif on the forward strand only.")

    g_data = p.add_argument_group("processed final dataset")
    g_data.add_argument("--final-dataset", default=None, help="CSV final-site table; bypasses all built-in raw-read processing.")
    g_data.add_argument("--candidate-model", choices=("auto", "motif", "all-bases", "observed"), default="auto", help="How omitted zero-count candidate sites are completed before classification.")
    g_data.add_argument("--classifier", default=None, help="Trusted local .py custom classifier exposing classify_gene(gene_id, site_rows).")
    g_data.add_argument("--no-essentiality", action="store_true", help="Do not call/colour gene essentiality.")
    g_data.add_argument("--no-read-histogram", action="store_true", help="Do not draw the blue per-site count histogram.")
    g_data.add_argument("--read-histogram-cap", type=float, default=None, metavar="PCT", help="Cap the read-histogram scale at this percentile (e.g. 95) of each contig's positive counts, so one hypersaturated site does not squash the rest; taller bars clip to full height and the scale reads '≥ N'.")
    g_data.add_argument("--csv-dir", default=None, help="Directory for normalised final-site and gene-call CSV outputs.")
    g_data.add_argument("--contig-alias", action="append", default=[], metavar="INPUT=GENBANK", help="Map an input CSV/WIG contig ID to the GenBank accession (repeatable).")

    g_name = p.add_argument_group("labelling")
    g_name.add_argument("--name", default=None, help="Override genome name/accession shown on plot.")
    g_name.add_argument("--no-orf-finder", action="store_true", help="Do not call de-novo ORFs when GenBank has no CDS.")
    g_name.add_argument("--no-legend", action="store_true", help="Hide legends.")
    g_name.add_argument("--no-read-legend", action="store_true", help="Hide the numeric read-count scale beside the histogram (shown by default).")
    g_name.add_argument("--small-title", action="store_true", help="Show a small label beside first row instead of large title.")

    g_tracks = p.add_argument_group("optional tracks")
    g_tracks.add_argument("--gc-content", action="store_true", help="Add GC-content deviation track.")
    g_tracks.add_argument("--gc-skew", action="store_true", help="Add GC-skew track.")
    g_tracks.add_argument("--trna", action="store_true", help="Draw tRNA/tmRNA/CRISPR features.")
    g_tracks.add_argument("--gc-window", type=int, default=None, help="GC window in bp; default is genome length/100.")

    g_out = p.add_argument_group("output size & style")
    g_out.add_argument("--paper", choices=sorted(PAPER_SIZES), default=None, help="Fit to a standard paper size.")
    g_out.add_argument("--portrait", action="store_true", help="Use portrait orientation for --paper.")
    g_out.add_argument("--fit-page", action="store_true", help="Force exact paper dimensions.")
    g_out.add_argument("--width", type=float, default=15.0, help="Figure width in inches.")
    g_out.add_argument("--track-height", type=float, default=0.93, help="Height per genome track in inches.")
    g_out.add_argument("--wrap-kb", type=float, default=20.0, help="Wrap genome rows about every N kb; 0 = one line.")
    g_out.add_argument("--rows", type=int, default=None, help="Force exact row count.")
    g_out.add_argument("--no-wrap", action="store_true", help="Draw the whole genome on one line.")
    g_out.add_argument("--dpi", type=int, default=300, help="PNG resolution.")
    g_out.add_argument("--transparent", action="store_true", help="Transparent PNG/SVG background.")


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    try:
        if argv and argv[0] == "process":
            return _run_process(build_process_parser().parse_args(argv[1:]))
        if argv and argv[0] == "plot":
            return _run_plot(build_plot_parser().parse_args(argv[1:]))
        return _run_plot(build_parser().parse_args(argv))
    except (DatasetError, PipelineError, ValueError, OSError, RuntimeError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1


def _run_plot(args: argparse.Namespace) -> int:
    reference = Path(args.input)
    _require_file(reference, "input genome")
    transposon = _resolve_transposon_or_raise(args.transposon)
    records = _load_records(reference, args.name, args.no_orf_finder)
    insertion_sites = _theoretical_sites(records, transposon, args)
    options = _plot_options(args, bool(insertion_sites))
    output = Path(args.output)

    if not args.final_dataset:
        out = render(records, insertion_sites, options, output, transposon_label=_transposon_label(transposon))
        print(f"Done. Wrote {out}")
        return 0

    default_contig = records[0].accession if len(records) == 1 else None
    sites = load_final_dataset(args.final_dataset, default_contig=default_contig)
    sites = _apply_contig_aliases(sites, _parse_contig_aliases(args.contig_alias))
    candidate_model = _resolve_candidate_model(args.candidate_model, transposon)
    if candidate_model != "observed":
        candidates = _candidate_sites(records, transposon, candidate_model, args.single_strand)
        sites = fill_missing_final_sites(sites, candidates)
    if candidate_model == "all-bases" and not args.no_essentiality:
        print(
            "warning: all-bases candidate model adapts a TA-site classifier; validate it for your transposon.",
            file=sys.stderr,
        )

    annotations, result = _annotate_and_classify(
        records, sites, classifier_path=args.classifier, skip_essentiality=args.no_essentiality
    )
    csv_dir = Path(args.csv_dir) if args.csv_dir else (output.parent if output.parent != Path("") else Path("."))
    site_csv, gene_csv = _write_analysis_csvs(
        sites, annotations, result, csv_dir, stem=output.stem
    )
    out = _render_final_dataset(
        records, insertion_sites, options, output, transposon, sites, result
    )
    print(f"Done. Wrote {out}\n      • {site_csv}")
    if gene_csv:
        print(f"      • {gene_csv}")
    return 0


def _run_process(args: argparse.Namespace) -> int:
    reference = Path(args.input)
    reads1 = Path(args.reads1)
    reads2 = Path(args.reads2) if args.reads2 else None
    _require_file(reference, "input genome")
    _require_file(reads1, "read 1")
    if reads2 is not None:
        _require_file(reads2, "read 2")
    if args.interactive:
        _interactive_process(args)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    replicon_ids = tuple(item.strip() for item in args.tpp_replicon_ids.split(",") if item.strip()) if args.tpp_replicon_ids else None
    tools = ToolPaths(
        fastqc=args.fastqc_bin,
        fastp=args.fastp_bin,
        tpp=(args.tpp_python, args.tpp_bin) if args.tpp_python else (args.tpp_bin,),
        bwa=args.bwa_bin,
        seqkit=args.seqkit_bin,
    )
    config = PipelineConfig(
        output_dir=output_dir,
        threads=args.threads,
        run_fastqc=not args.skip_fastqc,
        run_fastp=not args.skip_fastp,
        run_tpp=not args.skip_tpp,
        minimum_phred=args.quality_phred,
        minimum_read_length=args.min_read_length,
        adapter_sequence=args.adapter_sequence,
        adapter_sequence_r2=args.adapter_sequence_r2,
        tpp_primer=args.tpp_primer,
        tpp_mismatches=args.tpp_mismatches,
        tpp_mode=args.tpp_mode,
        tpp_replicon_ids=replicon_ids,
        minimum_mapped_reads=args.min_mapped_reads,
        include_unobserved_sites=not args.observed_sites_only,
        subsample_depth=args.subsample_depth,
        subsample_replicates=args.subsample_replicates,
        subsample_seed=args.subsample_seed,
    )
    print("[1/3] Processing reads (selected external stages only) ...")
    result = run_pipeline(
        reference,
        ReadPair(reads1, reads2, args.sample_id),
        config,
        tools=tools,
        runner=_logging_runner(output_dir / "pipeline.log"),
    )
    _write_manifest(output_dir / "processing_manifest.json", config, tools, result.commands, reference, reads1, reads2)

    counts = result.averaged_counts or result.counts
    if not counts:
        print("TPP was skipped: wrote reference/processing artifacts only; use `plot --final-dataset` to visualise externally processed counts.")
        return 0

    print("[2/3] Annotating final insertion sites and calling essentiality ...")
    records = _load_records(reference, None, False)
    sites = final_sites_from_counts(counts)
    annotations, classification = _annotate_and_classify(
        records, sites, classifier_path=args.classifier, skip_essentiality=args.no_essentiality
    )
    site_csv, gene_csv = _write_analysis_csvs(
        sites, annotations, classification, output_dir, stem="final"
    )
    print(f"      • final site table: {site_csv}")
    if gene_csv:
        print(f"      • gene calls: {gene_csv}")

    if args.no_plot:
        return 0
    print("[3/3] Rendering processed-data map ...")
    transposon = resolve_transposon("mariner" if args.tpp_mode == "himar1" else "tn5")
    if args.show_candidate_sites and transposon.has_preference:
        insertion_sites = {
            rec.accession: find_insertion_sites(rec.sequence, transposon.motif or "TA", both_strands=False)
            for rec in records
        }
    else:
        insertion_sites = {}
    plot_output = Path(args.plot_output)
    if not plot_output.is_absolute():
        plot_output = output_dir / plot_output
    options = PlotOptions(
        show_insertion_sites=bool(insertion_sites),
        show_insertion_density=False,
        show_read_histogram=not args.no_read_histogram,
        read_histogram_cap_percentile=_resolve_read_histogram_cap(args.read_histogram_cap),
    )
    out = _render_final_dataset(
        records, insertion_sites, options, plot_output, transposon, sites, classification
    )
    print(f"Done. Wrote {out}")
    return 0


def _load_records(reference: Path, name: str | None, no_orf_finder: bool) -> list[GenomeRecord]:
    print(f"Loading genome from {reference} ...")
    records = load_genome(reference, name_override=name, call_orfs_if_missing=not no_orf_finder)
    for record in records:
        source = "annotated GenBank" if record.annotation_source == "genbank" else "de-novo ORFs (pyrodigal-gv)"
        print(f"  • {record.accession}: {record.length:,} bp, {record.n_genes} genes [{source}]")
    return records


def _theoretical_sites(
    records: list[GenomeRecord], transposon: Transposon, args: argparse.Namespace
) -> dict[str, list[int]]:
    """Compute candidate insertion sites for the tick and/or density tracks.

    Returns the ``{accession: [positions]}`` mapping (empty when unavailable). Whether the
    red ticks and the density heat strip are actually *drawn* is decided independently in
    ``_plot_options`` — both are derived from this one mapping, so ``--no-insertion-sites``
    can hide the ticks while the density track stays on, and ``--no-insertion-density``
    the reverse. Only skip the computation when neither track will be drawn.
    """
    # A processed-data map is intentionally centred on measured blue read bars.
    # The dense red candidate-site tracks remain available, but are legacy
    # pre-visualisation output rather than a useful default at this stage.
    if args.final_dataset and not args.show_theoretical_sites:
        return {}
    if args.no_insertion_sites and args.no_insertion_density:
        return {}
    if not transposon.has_preference:
        if not args.no_insertion_sites:
            print(f"  • {transposon.name}: no finite motif preference; theoretical sites skipped.")
        return {}
    sites = {
        record.accession: find_insertion_sites(
            record.sequence, transposon.motif or "", both_strands=not args.single_strand
        )
        for record in records
    }
    for accession, positions in sites.items():
        print(f"  • {accession}: {len(positions):,} candidate {transposon.motif} sites")
    return sites


def _plot_options(args: argparse.Namespace, sites_available: bool) -> PlotOptions:
    return PlotOptions(
        fig_width=args.width,
        track_height=args.track_height,
        paper=args.paper,
        portrait=args.portrait,
        fit_page=args.fit_page,
        wrap_kb=0.0 if args.no_wrap else args.wrap_kb,
        force_rows=args.rows,
        dpi=args.dpi,
        transparent=args.transparent,
        show_insertion_sites=sites_available and not args.no_insertion_sites,
        show_insertion_density=sites_available and not args.no_insertion_density,
        show_read_histogram=not args.no_read_histogram,
        read_histogram_cap_percentile=_resolve_read_histogram_cap(args.read_histogram_cap),
        show_gc_content=args.gc_content,
        show_gc_skew=args.gc_skew,
        show_trna=args.trna,
        gc_window=args.gc_window,
        density_window=args.density_window,
        show_legend=not args.no_legend,
        show_read_legend=not args.no_read_legend,
        big_title=not args.small_title,
    )


def _resolve_read_histogram_cap(pct: float | None) -> float | None:
    if pct is None:
        return None
    if not 0 < pct <= 100:
        raise ValueError("--read-histogram-cap must be a percentile in (0, 100]")
    return pct


def _resolve_transposon_or_raise(spec: str) -> Transposon:
    try:
        return resolve_transposon(spec)
    except ValueError as exc:
        raise ValueError(f"invalid transposon: {exc}") from exc


def _transposon_label(transposon: Transposon) -> str:
    return transposon.motif if transposon.has_preference else transposon.name


def _resolve_candidate_model(model: str, transposon: Transposon) -> str:
    if model == "auto":
        return "motif" if transposon.has_preference else "all-bases"
    if model == "motif" and not transposon.has_preference:
        raise ValueError("--candidate-model motif requires a transposon with a finite IUPAC motif")
    return model


def _candidate_sites(
    records: Iterable[GenomeRecord], transposon: Transposon, model: str, single_strand: bool
) -> dict[str, list[int]]:
    if model == "observed":
        return {}
    candidates: dict[str, list[int]] = {}
    for record in records:
        if model == "all-bases":
            candidates[record.accession] = list(range(1, record.length + 1))
        else:
            candidates[record.accession] = find_insertion_sites(
                record.sequence, transposon.motif or "", both_strands=not single_strand
            )
    return candidates


def _parse_contig_aliases(specs: Iterable[str]) -> dict[str, str]:
    aliases: dict[str, str] = {}
    for spec in specs:
        if "=" not in spec:
            raise ValueError(f"invalid --contig-alias {spec!r}; expected INPUT=GENBANK")
        source, target = (part.strip() for part in spec.split("=", 1))
        if not source or not target:
            raise ValueError(f"invalid --contig-alias {spec!r}; neither side may be blank")
        if source in aliases and aliases[source] != target:
            raise ValueError(f"conflicting aliases supplied for contig {source!r}")
        aliases[source] = target
    return aliases


def _apply_contig_aliases(sites: Iterable[FinalSite], aliases: Mapping[str, str]) -> list[FinalSite]:
    return [replace(site, contig=aliases.get(site.contig, site.contig)) for site in sites]


def _annotate_and_classify(
    records: Iterable[GenomeRecord],
    sites: Iterable[FinalSite],
    *,
    classifier_path: str | None,
    skip_essentiality: bool,
) -> tuple[list[AnnotatedInsertionSite], ClassificationResult | None]:
    records = list(records)
    sites = list(sites)
    valid_contigs = {record.accession for record in records}
    unknown = sorted({site.contig for site in sites} - valid_contigs)
    if unknown:
        raise DatasetError(
            "final dataset contig(s) do not match GenBank accession(s): "
            f"{', '.join(unknown)}. Use --contig-alias if needed."
        )

    annotations: list[AnnotatedInsertionSite] = []
    for record in records:
        raw_sites = [
            InsertionSite(position=site.position, read_count=site.read_count, contig=site.contig)
            for site in sites
            if site.contig == record.accession
        ]
        annotations.extend(annotate_sites_with_genes(raw_sites, record.genes, contig=record.accession))

    if skip_essentiality:
        if classifier_path:
            raise ValueError("--classifier cannot be combined with --no-essentiality")
        return annotations, None
    classifier = load_classifier(classifier_path) if classifier_path else None
    return annotations, classify_genes(annotations, classifier=classifier)


def _write_analysis_csvs(
    sites: Iterable[FinalSite],
    annotations: Iterable[AnnotatedInsertionSite],
    classification: ClassificationResult | None,
    output_dir: Path,
    *,
    stem: str,
) -> tuple[Path, Path | None]:
    output_dir.mkdir(parents=True, exist_ok=True)
    assignments = {
        (site.contig, site.position): tuple((gene.gene_id, gene.strand) for gene in site.genes)
        for site in annotations
    }
    site_path = write_final_dataset(output_dir / f"{stem}_sites.csv", sites, gene_assignments=assignments)
    if classification is None:
        return site_path, None
    gene_path = output_dir / f"{stem}_gene_essentiality.csv"
    with gene_path.open("w", encoding="utf-8", newline="") as handle:
        fields = [
            "contig", "gene_id", "strand", "total_candidate_sites", "hit_sites",
            "saturation", "initial_call", "final_call", "read_count_median_threshold",
        ]
        import csv

        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for call in classification.calls:
            threshold = call.read_count_median_threshold
            writer.writerow(
                {
                    "contig": call.contig,
                    "gene_id": call.gene_id,
                    "strand": call.strand,
                    "total_candidate_sites": call.total_sites,
                    "hit_sites": call.hits,
                    "saturation": f"{call.saturation:g}",
                    "initial_call": call.initial_call or "",
                    "final_call": call.final_call or "",
                    "read_count_median_threshold": "" if threshold is None else f"{threshold:g}",
                }
            )
    return site_path, gene_path


def _render_final_dataset(
    records: list[GenomeRecord],
    insertion_sites: dict[str, list[int]],
    options: PlotOptions,
    output: Path,
    transposon: Transposon,
    sites: Iterable[FinalSite],
    classification: ClassificationResult | None,
) -> Path:
    calls = None
    if classification is not None:
        calls = {(call.contig, call.gene_id): call.final_call for call in classification.calls}
    return render(
        records,
        insertion_sites,
        options,
        output,
        transposon_label=_transposon_label(transposon),
        read_counts=group_counts_for_plotting(sites),
        gene_calls=calls,
    )


def _logging_runner(log_path: Path):
    log_path.parent.mkdir(parents=True, exist_ok=True)

    def run(command: tuple[str, ...]) -> None:
        command_text = shlex.join(command)
        print(f"  $ {command_text}")
        with log_path.open("a", encoding="utf-8") as log:
            log.write(f"\n$ {command_text}\n")
            log.flush()
            try:
                subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
            except FileNotFoundError as exc:
                raise PipelineError(
                    f"required external program was not found: {command[0]!r}; "
                    "install it or set its --*-bin path"
                ) from exc
            except subprocess.CalledProcessError as exc:
                raise PipelineError(
                    f"external command failed with exit code {exc.returncode}; see {log_path}"
                ) from exc

    return run


def _write_manifest(
    path: Path,
    config: PipelineConfig,
    tools: ToolPaths,
    commands: Iterable[tuple[str, ...]],
    reference: Path,
    reads1: Path,
    reads2: Path | None,
) -> None:
    payload = {
        "reference_genbank": str(reference),
        "reads1": str(reads1),
        "reads2": str(reads2) if reads2 else None,
        "parameters": _json_safe(asdict(config)),
        "tools": _json_safe(asdict(tools)),
        "commands": [list(command) for command in commands],
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _json_safe(value):
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_safe(item) for item in value]
    return value


def _interactive_process(args: argparse.Namespace) -> None:
    """Prompt only when explicitly requested, while retaining batch CLI flags."""
    args.skip_fastqc = not _ask_yes_no("Run FastQC initial quality assessment", not args.skip_fastqc)
    args.skip_fastp = not _ask_yes_no("Run fastp adapter trimming and quality filtering", not args.skip_fastp)
    if not args.skip_fastp:
        args.quality_phred = _ask_int("Minimum qualified Phred score", args.quality_phred, minimum=0)
        args.min_read_length = _ask_int("Minimum read length after trimming", args.min_read_length, minimum=1)
    args.skip_tpp = not _ask_yes_no("Run TRANSIT TPP + its integrated BWA mapping", not args.skip_tpp)
    if not args.skip_tpp:
        args.tpp_mismatches = _ask_int("TPP allowed transposon-prefix mismatches", args.tpp_mismatches, minimum=0)
        args.min_mapped_reads = _ask_float("Mapped-read count threshold (less than becomes zero)", args.min_mapped_reads, minimum=0)
        if _ask_yes_no("Depth-match with random SeqKit subsamples", args.subsample_depth is not None):
            args.subsample_depth = _ask_int("Reads per subsample", args.subsample_depth or 1, minimum=1)
            args.subsample_replicates = _ask_int("Number of subsamples to average", args.subsample_replicates, minimum=1)
        else:
            args.subsample_depth = None
            args.subsample_replicates = 1


def _ask_yes_no(question: str, default: bool) -> bool:
    prompt = "Y/n" if default else "y/N"
    while True:
        try:
            answer = input(f"{question} [{prompt}]: ").strip().lower()
        except EOFError as exc:
            raise PipelineError("interactive input ended unexpectedly") from exc
        if not answer:
            return default
        if answer in {"y", "yes"}:
            return True
        if answer in {"n", "no"}:
            return False
        print("Please answer y or n.")


def _ask_int(question: str, default: int, *, minimum: int) -> int:
    while True:
        try:
            value = int(input(f"{question} [{default}]: ").strip() or default)
        except (ValueError, EOFError):
            print(f"Enter an integer ≥ {minimum}.")
            continue
        if value >= minimum:
            return value
        print(f"Enter an integer ≥ {minimum}.")


def _ask_float(question: str, default: float, *, minimum: float) -> float:
    while True:
        try:
            value = float(input(f"{question} [{default:g}]: ").strip() or default)
        except (ValueError, EOFError):
            print(f"Enter a number ≥ {minimum:g}.")
            continue
        if value >= minimum:
            return value
        print(f"Enter a number ≥ {minimum:g}.")


def _require_file(path: Path, label: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"{label} not found: {path}")


if __name__ == "__main__":
    raise SystemExit(main())
