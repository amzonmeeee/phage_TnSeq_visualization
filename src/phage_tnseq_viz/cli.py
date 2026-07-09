"""Command-line interface for phage-tnseq-viz."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import __version__
from .genome import load_genome
from .plot import PAPER_SIZES, PlotOptions, render
from .transposon import PRESETS, find_insertion_sites, resolve_transposon


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="phage-tnseq-viz",
        description=(
            "Pre-visualization for phage Tn-Seq: draw a linear phage genome with "
            "gene arrows and transposon insertion sites (no sequencing data yet)."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("input", help="Input phage genome in GenBank (.gbk/.gb) format.")
    p.add_argument("-o", "--output", default="phage_map.png",
                   help="Output image path. Extension sets the format (.png or .svg).")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    g_tn = p.add_argument_group("transposon")
    g_tn.add_argument(
        "-t", "--transposon", default="mariner",
        help="Preset name (%s) or a custom IUPAC insertion motif (e.g. TA, NTAN)."
        % ", ".join(sorted(k for k in PRESETS)),
    )
    g_tn.add_argument("--no-insertion-sites", action="store_true",
                      help="Do not draw insertion-site ticks.")
    g_tn.add_argument("--insertion-density", action="store_true",
                      help="Add a sliding-window insertion-density heat track (sites/kb).")
    g_tn.add_argument("--density-window", type=int, default=None,
                      help="Window size (bp) for the density track. Default: length/100.")
    g_tn.add_argument("--single-strand", action="store_true",
                      help="Scan the forward strand only for the insertion motif.")

    g_name = p.add_argument_group("labelling")
    g_name.add_argument("--name", default=None,
                        help="Override the genome name/accession shown on the plot.")
    g_name.add_argument("--no-orf-finder", action="store_true",
                        help="Do not call ORFs de-novo when the GenBank has no CDS.")
    g_name.add_argument("--no-legend", action="store_true", help="Hide the legend.")
    g_name.add_argument("--small-title", action="store_true",
                        help="Show the genome name as a small label beside the first "
                             "row instead of a big heading on top.")

    g_tracks = p.add_argument_group("optional tracks")
    g_tracks.add_argument("--gc-content", action="store_true",
                          help="Add a GC-content deviation track.")
    g_tracks.add_argument("--gc-skew", action="store_true", help="Add a GC-skew track.")
    g_tracks.add_argument("--trna", action="store_true",
                          help="Draw tRNA/tmRNA/CRISPR features (from the GenBank).")
    g_tracks.add_argument("--gc-window", type=int, default=None,
                          help="Sliding-window size (bp) for GC tracks. Default: length/100.")

    g_out = p.add_argument_group("output size & style")
    g_out.add_argument("--paper", choices=sorted(PAPER_SIZES), default=None,
                       help="Fit to a standard paper size.")
    g_out.add_argument("--portrait", action="store_true",
                       help="Use portrait orientation for --paper.")
    g_out.add_argument("--fit-page", action="store_true",
                       help="Force the figure to the exact paper dimensions.")
    g_out.add_argument("--width", type=float, default=15.0,
                       help="Figure width in inches (ignored if --paper is set).")
    g_out.add_argument("--track-height", type=float, default=0.93,
                       help="Height per genome track (row) in inches.")
    g_out.add_argument("--wrap-kb", type=float, default=20.0,
                       help="Wrap the genome onto a new row about every N kb "
                            "(rows are evenly balanced). 0 = single line.")
    g_out.add_argument("--rows", type=int, default=None,
                       help="Force an exact number of rows (overrides --wrap-kb).")
    g_out.add_argument("--no-wrap", action="store_true",
                       help="Draw the whole genome on a single line.")
    g_out.add_argument("--dpi", type=int, default=300, help="Raster resolution for PNG.")
    g_out.add_argument("--transparent", action="store_true",
                       help="Transparent background (PNG/SVG).")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    in_path = Path(args.input)
    if not in_path.exists():
        print(f"error: input file not found: {in_path}", file=sys.stderr)
        return 2

    try:
        transposon = resolve_transposon(args.transposon)
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    print(f"[1/3] Loading genome from {in_path} ...")
    try:
        records = load_genome(
            in_path,
            name_override=args.name,
            call_orfs_if_missing=not args.no_orf_finder,
        )
    except Exception as exc:  # noqa: BLE001
        print(f"error: failed to load genome: {exc}", file=sys.stderr)
        return 1

    for rec in records:
        src = "annotated GenBank" if rec.annotation_source == "genbank" else "de-novo ORFs (pyrodigal-gv)"
        print(f"      • {rec.accession}: {rec.length:,} bp, {rec.n_genes} genes [{src}]")

    print(f"[2/3] Finding insertion sites for '{transposon.name}' ...")
    insertion_sites: dict[str, list[int]] = {}
    draw_insertions = not args.no_insertion_sites
    if draw_insertions and transposon.has_preference:
        for rec in records:
            sites = find_insertion_sites(
                rec.sequence, transposon.motif, both_strands=not args.single_strand
            )
            insertion_sites[rec.accession] = sites
            print(f"      • {rec.accession}: {len(sites):,} '{transposon.motif}' sites")
    elif draw_insertions and not transposon.has_preference:
        print(f"      • '{transposon.name}' has no strong sequence preference "
              f"(≈ random) — insertion ticks skipped.")
        draw_insertions = False

    options = PlotOptions(
        fig_width=args.width,
        track_height=args.track_height,
        paper=args.paper,
        portrait=args.portrait,
        fit_page=args.fit_page,
        wrap_kb=0.0 if args.no_wrap else args.wrap_kb,
        force_rows=args.rows,
        big_title=not args.small_title,
        dpi=args.dpi,
        transparent=args.transparent,
        show_insertion_sites=draw_insertions,
        show_insertion_density=args.insertion_density and draw_insertions,
        density_window=args.density_window,
        show_gc_content=args.gc_content,
        show_gc_skew=args.gc_skew,
        show_trna=args.trna,
        gc_window=args.gc_window,
        show_legend=not args.no_legend,
    )

    tn_label = transposon.motif if transposon.has_preference else transposon.name
    print(f"[3/3] Rendering -> {args.output} ...")
    try:
        out = render(records, insertion_sites, options, args.output,
                     transposon_label=tn_label or "")
    except Exception as exc:  # noqa: BLE001
        print(f"error: failed to render plot: {exc}", file=sys.stderr)
        return 1

    print(f"Done. Wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
