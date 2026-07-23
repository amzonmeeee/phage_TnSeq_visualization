"""Render a linear phage genome map with optional Tn-Seq data overlays.

Layout per genome/contig, top to bottom within one feature track:
  * gene arrows  – optional essentiality-coloured fill, border colour = PHROG
                   functional category (phold scheme); de-novo ORFs get a grey
                   dashed border.
  * optional blue read-count histogram immediately above each arrow lane.
  * insertion track – short red ticks at every transposon insertion site.
  * (optional) GC content, GC skew, tRNA/CRISPR sub-tracks.

Built on pyGenomeViz (matplotlib) so output can be PNG or SVG, transparent, and
sized to a custom figure or a standard paper size.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

# The package only writes figures; a non-interactive backend makes CLI rendering
# reliable on clusters and CI machines without a display server.
import matplotlib

matplotlib.use("Agg")

from pygenomeviz import GenomeViz

from . import colors
from .genome import GenomeRecord, gene_identifier
from .tracks import (
    WindowTrack,
    compute_gc_content,
    compute_gc_skew,
    compute_insertion_density,
)

# Paper sizes in inches, landscape (width, height).
PAPER_SIZES: dict[str, tuple[float, float]] = {
    "a5": (8.27, 5.83),
    "a4": (11.69, 8.27),
    "a3": (16.54, 11.69),
    "a2": (23.39, 16.54),
    "letter": (11.0, 8.5),
    "legal": (14.0, 8.5),
}

GC_POSITIVE = "#000000"   # GC content above genome mean
GC_NEGATIVE = "#808080"   # below mean
SKEW_POSITIVE = "#2ca02c" # green
SKEW_NEGATIVE = "#8900ff" # purple
TRNA_COLOR = "#111111"

# Gene arrow shaft height relative to the full track lane. 1.0 gives a "pentagon"
# arrow: a full-height body ending in a triangular point (home-plate shape), as
# opposed to the classic thin-tailed arrow. A gene clipped at a wrap boundary (the
# piece without the point) is drawn as a plain rectangle at exactly this height so
# it matches the neighbouring arrow bodies.
ARROW_SHAFT_RATIO = 1.0

# Arrow geometry and the data lanes matched to it. pyGenomeViz draws gene arrows
# across data y in [-ARROW_HALF, ARROW_HALF]. The insertion (red) tick lane, the
# density (heat) lane, and a maximal read-count bar are all sized so that a
# full-height element equals the *drawn* arrow height — the arrow, the red ticks and
# the heat strip end up the same thickness.
ARROW_HALF = 1.0
# Data-unit clearance kept beyond each arrow edge so the thick border stroke does not
# bleed into the neighbouring lane (the red ticks below / the histogram baseline above).
LANE_BORDER_PAD = 0.55
# The read-count histogram sits this far above the arrow top, and its tallest bar
# reaches READ_HIST_AMPLITUDE data units — 3 arrow heights (arrow height = 2 * ARROW_HALF).
READ_HIST_BASELINE_GAP = 0.55
READ_HIST_AMPLITUDE = 6.0
READ_HIST_LINEWIDTH = 0.9
# GC content/skew are kept at their previous absolute height while the arrow lane
# was halved, so their ratio (relative to the now-thin arrow lane) is doubled.
GC_TRACK_RATIO = 2.3

# Gene-arrow lane height (pyGenomeViz feature_track_ratio). Halved from the 0.25
# default so the arrow / insertion / density lanes read as slim bars.
FEATURE_TRACK_RATIO = 0.125
# Gene arrow border line width.
GENE_EDGE_WIDTH = 4.95
# Target drawn height of the gene-arrow lane, as a share of one track's height. The
# lane's natural height is whatever the margins leave over, so it swings with the number
# of sub-tracks: without a floor the border stroke eats half the arrow (fat pills) on a
# figure with no sub-tracks, and a sixth of it when every track is on. Pinning the lane
# keeps the arrows — and the share of them taken by the border — the same either way.
MIN_GENE_LANE_RATIO = 0.475
# Bottom legend band, in figure fractions: the gap between cells, the width of the
# colorbar's own cell, and the width of the colour strip drawn inside that cell.
BAND_GAP = 0.018
BAND_LEFT, BAND_RIGHT = 0.04, 0.97
# Cells per band row. A leftover odd cell spans the whole row on its own, so that every
# row starts and ends on the same two edges.
BAND_COLS = 2
# Vertical gap between band rows, and the cell height used when the colorbar stands
# alone and so has no legend box to take its height from. Both in inches.
BAND_ROW_GAP_IN = 0.18
CBAR_ONLY_HEIGHT_IN = 1.25
# Legend text size, held fixed: cells take the widest group's natural width rather
# than shrinking the type to fit, and the row spills past the figure edges if need be.
LEGEND_FONTSIZE = 9.0
TITLE_FONTSIZE_RATIO = 9.5 / 8.0
FRAME_EDGE_COLOR = "#999999"


@dataclass
class PlotOptions:
    # figure sizing
    fig_width: float = 15.0
    track_height: float = 0.93
    paper: str | None = None
    portrait: bool = False
    fit_page: bool = False
    # line wrapping (multi-row layout)
    wrap_kb: float = 20.0             # start a new row roughly every N kb (0 = off)
    force_rows: int | None = None     # force an exact number of rows
    # output
    dpi: int = 300
    transparent: bool = False
    # optional tracks
    show_insertion_sites: bool = True
    show_insertion_density: bool = True
    show_read_histogram: bool = True
    # Cap the read-histogram scale at this percentile of a contig's positive
    # counts instead of its single maximum. Real Tn-Seq counts are heavy-tailed,
    # so one hypersaturated site otherwise squashes every other bar to a sliver.
    # Bars above the cap are clipped to full height and the scale reads "≥ N".
    # None keeps the previous max-based scaling.
    read_histogram_cap_percentile: float | None = None
    show_gc_content: bool = False
    show_gc_skew: bool = False
    show_trna: bool = False
    gc_window: int | None = None
    density_window: int | None = None
    density_cmap: str = "inferno"
    # styling
    insertion_color: str = "#e50000"
    read_histogram_color: str = "#1677ff"
    gene_fill: str = "#000000"       # used when no essentiality calls are supplied
    show_legend: bool = True
    show_read_legend: bool = True    # numeric read-count scale beside the histogram
    big_title: bool = True           # genome name as a big heading on top
                                     # (False = small label beside the first row)


def _gene_edge(gene) -> tuple[str, str]:
    """Return (edgecolor, linestyle) for a gene arrow."""
    if not gene.annotated:
        return colors.UNANNOTATED_EDGE, "--"
    if gene.is_vfdb_amr:
        return colors.VFDB_AMR_COLOR, "-"
    return colors.category_color(gene.function), "-"


def _gene_fill(accession: str, gene, gene_calls, default: str) -> str:
    """Resolve one arrow fill without changing legacy no-data rendering."""
    if gene_calls is None:
        return default
    return colors.essentiality_color(
        gene_calls.get((accession, gene_identifier(gene, contig=accession)))
    )


def _row_windows(length: int, wrap_kb: float, force_rows: int | None) -> list[tuple[int, int]]:
    """Split a genome of ``length`` bp into contiguous, balanced row windows.

    Returns ``[(start, end), ...]`` (0-based, half-open) covering ``[0, length)``.
    A single-row genome yields ``[(0, length)]`` so short genomes stay on one line.
    ``force_rows`` overrides the automatic wrap-width choice; ``wrap_kb <= 0``
    disables wrapping. Rows are made as even as possible (rather than fixed-width
    with a tiny trailing remainder) so every row is close to full width.
    """
    if length <= 0:
        return [(0, max(length, 0))]
    if force_rows and force_rows > 0:
        n = force_rows
    else:
        wrap_bp = int(wrap_kb * 1000) if wrap_kb and wrap_kb > 0 else 0
        n = 1 if (wrap_bp <= 0 or length <= wrap_bp) else math.ceil(length / wrap_bp)
    n = max(1, min(n, length))
    width = math.ceil(length / n)
    windows: list[tuple[int, int]] = []
    start = 0
    while start < length and len(windows) < n:
        end = min(start + width, length)
        windows.append((start, end))
        start = end
    return windows or [(0, length)]


def _fmt_kb_range(r0: int, r1: int) -> str:
    return f"{r0 / 1000:g}–{r1 / 1000:g} kb"


def _row_label(rec: GenomeRecord, idx: int, r0: int, r1: int, n_rows: int,
               big_title: bool) -> str:
    # With a big top title the genome identity lives in the heading, so the rows
    # only carry their kb range (nothing at all when there is a single row).
    if big_title:
        return "" if n_rows == 1 else _fmt_kb_range(r0, r1)
    ident = f"{rec.name}\n{rec.accession} | {rec.length:,} bp"
    if n_rows == 1:
        return ident
    if idx == 0:
        return f"{ident}\n{_fmt_kb_range(r0, r1)}"
    return _fmt_kb_range(r0, r1)


def _draw_clipped_arrow(ax, gene, r0: int, row_width: int, head_length: float,
                        edge: str, linestyle: str, fill: str) -> None:
    """Draw a gene's *whole* arrow shape, clipped to one wrapped row.

    A boundary-straddling gene is drawn at its true full extent (arrowhead only at
    the real 3' end) and then clipped to this row's window. The clip spans x only:
    it cuts the wrap edge open — so the gene reads as one complete arrow sliced by
    the wrap line, continuing on the neighbouring row — while leaving y unbounded so
    the top/bottom borders keep their full line width. (Clipping to the arrow's own
    y-extent would halve them, since the stroke is centred on the polygon edge, and
    pyGenomeViz draws the unclipped neighbours with no clipping at all.)
    """
    from matplotlib.patches import Polygon
    from matplotlib.transforms import Bbox, TransformedBbox

    g0, g1 = gene.start - 1, gene.end
    xs, xe = g0 - r0, g1 - r0            # arrow spans axis-x [xs, xe] (may exceed row)
    y_lo, y_hi = ax.get_ylim()
    ym = (y_lo + y_hi) / 2.0
    hl = min(head_length, xe - xs)      # same head length pyGenomeViz uses
    if gene.strand >= 0:                 # head at the 3' (right) end
        verts = [(xs, y_lo), (xs, y_hi), (xe - hl, y_hi), (xe, ym), (xe - hl, y_lo)]
    else:                                # head at the 3' (left) end
        verts = [(xe, y_lo), (xe, y_hi), (xs + hl, y_hi), (xs, ym), (xs + hl, y_lo)]
    poly = Polygon(verts, closed=True, facecolor=fill, edgecolor=edge,
                   linewidth=GENE_EDGE_WIDTH, linestyle=linestyle,
                   joinstyle="miter", zorder=3)
    ax.add_patch(poly)
    pad = (y_hi - y_lo) * 10.0  # >> the stroke's half-width, so y never clips
    poly.set_clip_box(TransformedBbox(
        Bbox.from_extents(0, y_lo - pad, row_width, y_hi + pad), ax.transData))
    poly.set_clip_on(True)


def render(
    records: list[GenomeRecord],
    insertion_sites: dict[str, list[int]],
    options: PlotOptions,
    out_path: str | Path,
    transposon_label: str = "",
    *,
    read_counts: dict[str, dict[int, float]] | None = None,
    gene_calls: dict[tuple[str, str], str | None] | None = None,
) -> Path:
    """Render a map and optional final-dataset overlays.

    ``read_counts`` maps ``accession -> {1-based insertion position: count}`` and
    is rendered as a blue, contig-normalised histogram above the gene arrows.
    ``gene_calls`` maps ``(accession, gene identifier)`` to calls produced by the
    built-in or a user classifier.  The identifier is the GenBank ``locus_tag``
    (or ``ID``/``label`` when loaded) and falls back to ``start-end:strand``.
    """
    out_path = Path(out_path)
    read_counts = read_counts or {}

    fig_width = options.fig_width
    if options.paper:
        pw, ph = PAPER_SIZES[options.paper]
        if options.portrait:
            pw, ph = ph, pw
        fig_width = pw

    gv = GenomeViz(
        fig_width=fig_width,
        fig_track_height=options.track_height,
        track_align_type="left",
        feature_track_ratio=FEATURE_TRACK_RATIO,
        # Wrapped rows are stacked feature tracks; pyGenomeViz would otherwise
        # insert a full-height "link" band between them (meant for genome-to-genome
        # alignment ribbons). Shrink it to a thin inter-row gap.
        link_track_ratio=0.10,
    )

    # Pre-compute the whole-genome window signals (pure numeric, no axes needed).
    # Each is stored with a genome-wide max so every wrapped row shares one scale.
    density_signal: dict[str, tuple[WindowTrack, float]] = {}
    gcc_signal: dict[str, tuple[WindowTrack, float]] = {}
    gcs_signal: dict[str, tuple[WindowTrack, float]] = {}
    read_count_max: dict[str, float] = {}
    for rec in records:
        if options.show_insertion_density and insertion_sites.get(rec.accession):
            win = options.density_window or max(200, rec.length // 100)
            wt = compute_insertion_density(
                insertion_sites[rec.accession], rec.length, win, max(1, win // 2)
            )
            density_signal[rec.accession] = (wt, max(wt.values, default=0.0))
        if options.show_gc_content:
            win = options.gc_window or _default_window(rec.length)
            wt = compute_gc_content(rec.sequence, win, _gc_step(win))
            gcc_signal[rec.accession] = (wt, max((abs(v) for v in wt.values), default=1.0) or 1.0)
        if options.show_gc_skew:
            win = options.gc_window or _default_window(rec.length)
            wt = compute_gc_skew(rec.sequence, win, _gc_step(win))
            gcs_signal[rec.accession] = (wt, max((abs(v) for v in wt.values), default=1.0) or 1.0)
        positive_counts = [count for count in read_counts.get(rec.accession, {}).values() if count > 0]
        if positive_counts:
            if options.read_histogram_cap_percentile is not None:
                # Never let the cap fall below the smallest positive count, or an
                # all-equal contig would divide by zero / clip everything.
                read_count_max[rec.accession] = max(
                    _percentile(positive_counts, options.read_histogram_cap_percentile),
                    min(positive_counts),
                )
            else:
                read_count_max[rec.accession] = max(positive_counts)
    # The scale reads "≥ N" only when the cap actually clips at least one site.
    read_hist_capped = any(
        count > read_count_max.get(accession, 0.0)
        for accession, counts in read_counts.items()
        for count in counts.values()
    )

    # Vertical layout of the gene-arrow axis (data units). The arrow spans
    # [-ARROW_HALF, ARROW_HALF]; a border-clearance pad is kept on each side, and the
    # read histogram (when present) claims headroom above. `arrow_fraction` is the share
    # of the axis the arrow itself takes, which the red/heat lanes are matched to so all
    # three read the same thickness.
    has_read_histogram = bool(options.show_read_histogram and read_count_max)
    arrow_y_lo = -(ARROW_HALF + LANE_BORDER_PAD)
    if has_read_histogram:
        read_hist_base = ARROW_HALF + READ_HIST_BASELINE_GAP
        arrow_y_hi = read_hist_base + READ_HIST_AMPLITUDE
    else:
        read_hist_base = None
        arrow_y_hi = ARROW_HALF + LANE_BORDER_PAD
    arrow_fraction = (2.0 * ARROW_HALF) / (arrow_y_hi - arrow_y_lo)

    # Keep references (keyed by genome + row index) so we can draw on the
    # sub-track axes after plotfig(), and remember each row's coordinate window.
    insertion_subtracks: dict[tuple[str, int], object] = {}
    density_subtracks: dict[tuple[str, int], object] = {}
    gc_content_subtracks: dict[tuple[str, int], object] = {}
    gc_skew_subtracks: dict[tuple[str, int], object] = {}
    row_window: dict[tuple[str, int], tuple[int, int]] = {}
    row_tracks: dict[tuple[str, int], object] = {}
    # Genes straddling a wrap boundary — drawn by hand after plotfig() as one whole
    # arrow clipped to each row (pyGenomeViz rejects out-of-range coordinates).
    straddling: list[tuple[tuple[str, int], object, str, str, str]] = []
    present_categories: set[str] = set()
    present_essentiality_calls: set[str | None] = set()
    any_denovo = False
    any_vfdb = False

    for rec in records:
        rows = _row_windows(rec.length, options.wrap_kb, options.force_rows)
        n_rows = len(rows)
        for idx, (r0, r1) in enumerate(rows):
            key = (rec.accession, idx)
            row_window[key] = (r0, r1)
            # Left-aligned segment spanning this row's true coordinate window; all
            # rows share one bp-per-inch scale so a shorter last row ends early.
            track = gv.add_feature_track(
                _row_label(rec, idx, r0, r1, n_rows, options.big_title),
                (r0, r1), labelsize=11,
            )
            row_tracks[key] = track
            seg = track.get_segment()

            for gene in rec.genes:
                g0, g1 = gene.start - 1, gene.end  # 0-based half-open
                if g1 <= r0 or g0 >= r1:
                    continue
                edge, ls = _gene_edge(gene)
                fill = _gene_fill(rec.accession, gene, gene_calls, options.gene_fill)
                if gene_calls is not None:
                    present_essentiality_calls.add(
                        gene_calls.get((rec.accession, gene_identifier(gene, contig=rec.accession)))
                    )
                if gene.annotated and not gene.is_vfdb_amr:
                    present_categories.add(colors.category_key(gene.function))
                any_denovo = any_denovo or (not gene.annotated)
                any_vfdb = any_vfdb or gene.is_vfdb_amr
                if g0 >= r0 and g1 <= r1:
                    # Fully inside this row: let pyGenomeViz draw the arrow.
                    seg.add_feature(
                        g0, g1, gene.strand,
                        plotstyle="bigarrow",
                        arrow_shaft_ratio=ARROW_SHAFT_RATIO,
                        fc=fill,
                        ec=edge,
                        lw=GENE_EDGE_WIDTH,
                        linestyle=ls,
                    )
                else:
                    # Straddles a wrap boundary: draw the whole arrow ourselves and
                    # clip it to this row so it reads as one arrow sliced by the wrap.
                    straddling.append((key, gene, edge, ls, fill))

            # optional tRNA / CRISPR as thin boxes on the main lane
            if options.show_trna:
                for f in rec.noncds:
                    f0, f1 = f.start - 1, f.end
                    if f1 <= r0 or f0 >= r1:
                        continue
                    seg.add_feature(
                        max(f0, r0), min(f1, r1), f.strand,
                        plotstyle="bigrbox", fc=TRNA_COLOR, ec=TRNA_COLOR, lw=0.5,
                    )

            if options.show_insertion_sites and insertion_sites.get(rec.accession):
                insertion_subtracks[key] = track.add_subtrack(
                    f"ins:{rec.accession}:{idx}", ratio=arrow_fraction, ylim=(0, 1)
                )
            if rec.accession in density_signal:
                density_subtracks[key] = track.add_subtrack(
                    f"dens:{rec.accession}:{idx}", ratio=arrow_fraction, ylim=(0, 1)
                )
            if options.show_gc_content:
                gc_content_subtracks[key] = track.add_subtrack(
                    f"gcc:{rec.accession}:{idx}", ratio=GC_TRACK_RATIO, ylim=(-1, 1)
                )
            if options.show_gc_skew:
                gc_skew_subtracks[key] = track.add_subtrack(
                    f"gcs:{rec.accession}:{idx}", ratio=GC_TRACK_RATIO, ylim=(-1, 1)
                )

    gv.set_scale_xticks(labelsize=10, ymargin=0.3)
    fig = gv.plotfig()

    # ---- boundary-straddling genes: one whole arrow, clipped to each row ----
    if straddling:
        max_size = max(t.ax.get_xlim()[1] for t in row_tracks.values())
        head_length = max_size * 0.015  # pyGenomeViz's bigarrow head length
        for key, gene, edge, ls, fill in straddling:
            r0, r1 = row_window[key]
            _draw_clipped_arrow(
                row_tracks[key].ax, gene, r0, r1 - r0, head_length,
                edge, ls, fill,
            )

    # ---- draw signals on sub-track axes (only available after plotfig) ----
    # Set the shared vertical range on every gene-arrow axis: pad below/above the arrow
    # so its thick border clears the neighbouring lanes, and (with sequencing data) make
    # room above for the read histogram. All rows share one range so arrows line up.
    for track in row_tracks.values():
        track.ax.set_ylim(arrow_y_lo, arrow_y_hi)

    # Keep the read histogram in the same main axis as its arrows: the headroom above
    # the arrow top gives a lane *above* the arrows even for wrapped genomes, rather than
    # a pyGenomeViz subtrack (which is always placed below). The tallest bar reaches one
    # full arrow height so the busiest sites are legible instead of a thin fuzz.
    if has_read_histogram:
        for key, track in row_tracks.items():
            accession, _ = key
            r0, r1 = row_window[key]
            values = [
                (position, count)
                for position, count in read_counts.get(accession, {}).items()
                if r0 < position <= r1 and count > 0
            ]
            if not values:
                continue
            vmax = read_count_max[accession]
            xs = [track.transform_coord(position) for position, _ in values]
            # min(..., 1.0) clips a count above the percentile cap to full height
            # so it cannot overshoot into the arrow lane above.
            tops = [read_hist_base + READ_HIST_AMPLITUDE * min(count / vmax, 1.0) for _, count in values]
            track.ax.vlines(
                xs, read_hist_base, tops, color=options.read_histogram_color,
                lw=READ_HIST_LINEWIDTH, zorder=5,
            )
        # A small vertical scale on the first row of each contig telling the reader how
        # many reads the tallest bar stands for (the bars are otherwise unitless). This is
        # the read-count "legend"; toggle it off with show_read_legend.
        if options.show_read_legend:
            for accession, vmax in read_count_max.items():
                scale_track = row_tracks.get((accession, 0))
                if scale_track is not None:
                    _draw_read_scale(
                        scale_track.ax, read_hist_base, arrow_y_hi, vmax,
                        options.read_histogram_color, capped=read_hist_capped,
                    )

    density_mappable = None
    for key, sub in insertion_subtracks.items():
        acc, _ = key
        r0, r1 = row_window[key]
        xs = [sub.transform_coord(p) for p in insertion_sites.get(acc, []) if r0 < p <= r1]
        if xs:
            sub.ax.vlines(xs, 0.0, 1.0, color=options.insertion_color, lw=0.4)

    for key, sub in density_subtracks.items():
        acc, _ = key
        wt, vmax = density_signal[acc]
        m = _draw_density_row(sub, wt, vmax, row_window[key], cmap=options.density_cmap)
        density_mappable = density_mappable or m

    for key, sub in gc_content_subtracks.items():
        acc, _ = key
        wt, vmax = gcc_signal[acc]
        _draw_window_row(sub, wt, vmax, row_window[key], GC_POSITIVE, GC_NEGATIVE)

    for key, sub in gc_skew_subtracks.items():
        acc, _ = key
        wt, vmax = gcs_signal[acc]
        _draw_window_row(sub, wt, vmax, row_window[key], SKEW_POSITIVE, SKEW_NEGATIVE)

    # ---- title, legends and colorbar, all kept off to the top/bottom margins ----
    _decorate(
        fig, records, options,
        present_categories=present_categories,
        include_vfdb=any_vfdb, include_denovo=any_denovo,
        transposon_label=transposon_label,
        insertion_site_count=sum(len(sites) for sites in insertion_sites.values()),
        essentiality_calls=present_essentiality_calls if gene_calls is not None else None,
        has_read_histogram=has_read_histogram,
        density_mappable=density_mappable,
        arrow_fraction=arrow_fraction,
    )

    if options.paper and options.fit_page:
        pw, ph = PAPER_SIZES[options.paper]
        if options.portrait:
            pw, ph = ph, pw
        fig.set_size_inches(pw, ph)

    fig.savefig(
        out_path,
        dpi=options.dpi,
        transparent=options.transparent,
        bbox_inches="tight",
    )
    return out_path


def _percentile(values: list[float], pct: float) -> float:
    """Linear-interpolated percentile (``pct`` in [0, 100]) of ``values``."""
    if not values:
        return 0.0
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    rank = (max(0.0, min(100.0, pct)) / 100.0) * (len(ordered) - 1)
    low = math.floor(rank)
    high = math.ceil(rank)
    if low == high:
        return ordered[low]
    return ordered[low] + (rank - low) * (ordered[high] - ordered[low])


def _default_window(length: int) -> int:
    return max(100, length // 100)


def _gc_step(window: int) -> int:
    """Sliding step for GC tracks — a small fraction of the window so the curve is
    finely sampled and smooth (phold-style) rather than coarsely polygonal."""
    return max(1, window // 10)


def _draw_read_scale(ax, base: float, top: float, vmax: float, color: str,
                     *, capped: bool = False) -> None:
    """Draw a small "0 – N reads" vertical scale just right of a histogram row.

    The read bars are drawn normalised (tallest = the contig's peak, or the
    percentile cap when one is set); this bar tells the reader what that top is in
    absolute reads, so the histogram carries a real unit.  When ``capped`` the top
    label reads ``≥ N`` because taller sites were clipped to full height.  Placed
    in blended coords (x in axes fraction, y in data) so it hugs the right spine
    (an always-clear margin — the left one carries pyGenomeViz's row label)
    regardless of the genome's bp scale.
    """
    trans = ax.get_yaxis_transform()
    x = 1.006
    tick = 0.004
    ax.plot([x, x], [base, top], transform=trans, color=color, lw=1.3,
            clip_on=False, zorder=6)
    for y in (base, top):
        ax.plot([x - tick, x + tick], [y, y], transform=trans, color=color, lw=1.3,
                clip_on=False, zorder=6)
    magnitude = f"{int(round(vmax)):,}" if vmax >= 1 else f"{vmax:g}"
    label = f"≥{magnitude}" if capped else magnitude
    ax.text(x + 2 * tick, top, label, transform=trans, ha="left", va="center",
            fontsize=7.5, color=color, clip_on=False)
    ax.text(x + 2 * tick, base, "0", transform=trans, ha="left", va="center",
            fontsize=7.5, color=color, clip_on=False)
    ax.text(x + 8 * tick, (base + top) / 2, "reads", transform=trans, ha="center",
            va="center", rotation=90, fontsize=7.5, color=color, clip_on=False)


def _draw_window_row(sub, wt: WindowTrack, vmax: float, window, pos_c, neg_c) -> None:
    """Draw one row's slice of a GC content/skew signal, scaled by a shared vmax."""
    r0, r1 = window
    pairs = [(p, v) for p, v in zip(wt.positions, wt.values) if r0 <= p <= r1]
    if not pairs:
        return
    xs = [sub.transform_coord(min(max(p, r0), r1)) for p, _ in pairs]
    ys = [v / vmax for _, v in pairs]  # scale into (-1, 1) with a genome-wide vmax
    sub.ax.fill_between(xs, ys, 0, where=[v >= 0 for v in ys], color=pos_c, lw=0)
    sub.ax.fill_between(xs, ys, 0, where=[v < 0 for v in ys], color=neg_c, lw=0)
    sub.ax.axhline(0, color="#00000055", lw=0.4)


def _draw_density_row(sub, wt: WindowTrack, vmax: float, window, *, cmap):
    """Draw one row's slice of the insertion-density heat strip.

    Uses a shared ``vmax`` so colours are comparable across wrapped rows; returns
    the mappable so a single colorbar can annotate the whole figure.
    """
    import numpy as np

    r0, r1 = window
    vals = [v for p, v in zip(wt.positions, wt.values) if r0 <= p <= r1]
    if not vals:
        return None
    data = np.array(vals).reshape(1, -1)
    im = sub.ax.imshow(
        data,
        extent=(sub.transform_coord(r0), sub.transform_coord(r1), 0, 1),
        aspect="auto", cmap=cmap, origin="lower", interpolation="nearest",
        vmin=0.0, vmax=vmax or None,
    )
    return im


def _gene_lane_height_frac(fig) -> float:
    """Figure fraction taken by the gene-arrow lane, before any margins are reserved.

    Identified by the arrow stroke width, which nothing else on the figure uses.
    """
    for ax in fig.axes:
        for coll in ax.collections:
            lw = coll.get_linewidth()
            lw = lw[0] if hasattr(lw, "__len__") and len(lw) else lw
            if abs(float(lw) - GENE_EDGE_WIDTH) < 1e-6:
                return ax.get_position().height
    return fig.axes[0].get_position().height if fig.axes else 1.0


def _compress_axes(fig, top_in: float, bottom_in: float) -> tuple[float, float]:
    """Squeeze every existing axis into the vertical band [bottom, 1-top].

    Frees ``top_in`` / ``bottom_in`` inches at the figure edges for a title and a
    legend/colorbar strip, keeping the tracks' relative proportions intact.
    Returns the freed margins as figure fractions ``(bottom, top)``.
    """
    height = fig.get_size_inches()[1]
    t = top_in / height
    b = bottom_in / height
    span = max(1e-3, 1.0 - t - b)
    # Read every current position first, then apply — some axes are position-linked
    # (the scale ruler is a twin of the lowest track), so reading a position after
    # a sibling was moved would compound the shift.
    targets = [
        (ax, ax.get_position()) for ax in fig.axes
    ]
    for ax, p in targets:
        ax.set_position([p.x0, b + p.y0 * span, p.width, p.height * span])
    return b, t


def _pad_to_rows(handles: list, labels: list[str], ncol: int, rows: int) -> tuple[list, list]:
    """Pad a legend group with invisible entries so it lays out exactly ``rows`` rows.

    Every group then has the same number of rows and therefore the same framed height,
    which is what lets the boxes align into a grid. Matplotlib fills a legend column by
    column, so padding to ``rows * ncol`` entries keeps the columns balanced too.
    """
    from matplotlib.lines import Line2D

    handles, labels = list(handles), list(labels)
    while len(handles) < rows * ncol:
        handles.append(Line2D([], [], linestyle="none"))
        labels.append("")
    return handles, labels


def _draw_colorbar_cell(fig, mappable, x: float, top: float, width: float, height: float,
                        fs: float) -> None:
    """Draw the density colorbar as a framed cell matching the legend boxes.

    The frame, title and padding mimic a legend box so the colorbar reads as the last
    cell of the same grid rather than a loose annotation beside it.
    """
    from matplotlib.patches import FancyBboxPatch
    from matplotlib.ticker import MaxNLocator

    title_fs = fs * TITLE_FONTSIZE_RATIO

    fig.add_artist(FancyBboxPatch(
        (x, top - height), width, height,
        boxstyle="round,pad=0,rounding_size=0.008", transform=fig.transFigure,
        facecolor="white", edgecolor=FRAME_EDGE_COLOR, linewidth=0.8, zorder=4,
    ))
    fig.text(x + width / 2.0, top - 0.10 * height, "Insertion density (sites / kb)",
             ha="center", va="top", fontsize=title_fs, zorder=5)

    # Lay the bar out as shares of the cell, not fixed inches: the cell's height tracks
    # the legend boxes, so absolute padding would swallow the bar on short figures. The
    # bar sits above its tick labels, which need the lower third of the cell.
    cax = fig.add_axes([x + 0.10 * width, top - height + 0.38 * height,
                        0.80 * width, 0.20 * height], zorder=5)
    cbar = fig.colorbar(mappable, cax=cax, orientation="horizontal")
    cbar.outline.set_edgecolor(FRAME_EDGE_COLOR)
    cbar.locator = MaxNLocator(nbins=6, integer=True)
    cbar.update_ticks()
    cbar.ax.tick_params(labelsize=fs * 0.94, length=2.5, pad=2)


def _decorate(fig, records, options, *, present_categories, include_vfdb,
              include_denovo, transposon_label, insertion_site_count,
              essentiality_calls, has_read_histogram, density_mappable,
              arrow_fraction=1.0):
    """Add the top title plus framed legend boxes and a density colorbar below.

    All decorations live in reserved top/bottom margins (never beside the plot),
    and the gene-function legend is kept separate from the transposon/track one.
    """
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    rec = records[0]
    big = options.big_title

    # Gene-function (PHROG category) legend — border colour = category.
    gene_h: list = []
    gene_l: list[str] = []
    for label, color in colors.legend_entries(present_categories, include_vfdb=include_vfdb):
        gene_h.append(Patch(facecolor="white", edgecolor=color, linewidth=1.6))
        gene_l.append(label)
    if include_denovo:
        gene_h.append(Patch(facecolor="white", edgecolor=colors.UNANNOTATED_EDGE,
                            linewidth=1.6, linestyle="--"))
        gene_l.append("De-novo ORF (unannotated)")

    # Essentiality legend — arrow fill colour.  Calls are kept independent from
    # the PHROG border legend so both visual encodings remain readable.
    ess_h: list = []
    ess_l: list[str] = []
    if essentiality_calls is not None:
        for label, color in colors.essentiality_legend_entries(essentiality_calls):
            ess_h.append(Patch(facecolor=color, edgecolor="#555555", linewidth=0.8))
            ess_l.append(label)

    # Transposon legend (its own box).
    tn_h: list = []
    tn_l: list[str] = []
    if options.show_insertion_sites:
        motif = f" ({transposon_label})" if transposon_label else ""
        lbl = f"Insertion sites{motif} = {insertion_site_count:,}"
        tn_h.append(Line2D([0], [0], color=options.insertion_color, lw=1.6))
        tn_l.append(lbl)
    if has_read_histogram:
        tn_h.append(Line2D([0], [0], color=options.read_histogram_color, lw=2.0))
        tn_l.append("Read count per insertion site (blue histogram)")
    if essentiality_calls is None:
        tn_h.append(Patch(facecolor="black", edgecolor="#000000", linewidth=1.0))
        tn_l.append("Gene fill = essentiality (no data yet)")

    # With sequencing data, keep the count-bar/transposon key with the
    # essentiality palette.  It prevents a third, mostly-empty legend row while
    # still keeping PHROG borders in their own independent box.
    essentiality_title = "Tn-Seq essentiality"
    if ess_h and tn_h:
        ess_h.extend(tn_h)
        ess_l.extend(tn_l)
        tn_h, tn_l = [], []
        essentiality_title = "Tn-Seq data & essentiality"

    # Sequence-track legend (its own box).
    seq_h: list = []
    seq_l: list[str] = []
    if options.show_gc_content:
        seq_h.append(Patch(facecolor=GC_POSITIVE)); seq_l.append("GC content > genome mean")
        seq_h.append(Patch(facecolor=GC_NEGATIVE)); seq_l.append("GC content < genome mean")
    if options.show_gc_skew:
        seq_h.append(Patch(facecolor=SKEW_POSITIVE)); seq_l.append("GC skew +  (G > C)")
        seq_h.append(Patch(facecolor=SKEW_NEGATIVE)); seq_l.append("GC skew −  (C > G)")
    if options.show_trna:
        seq_h.append(Patch(facecolor=TRNA_COLOR)); seq_l.append("tRNA / tmRNA / CRISPR")

    # Each non-empty group becomes its own framed box; the boxes and the colorbar are
    # then laid out as equal cells of a grid, wrapping to a second row if one would
    # not fit across the figure.
    boxes: list[tuple[str, list, list, int]] = []
    if gene_h:
        boxes.append(("Gene function — PHROG category", gene_h, gene_l,
                      2 if len(gene_h) > 4 else 1))
    if ess_h:
        boxes.append((essentiality_title, ess_h, ess_l,
                      2 if len(ess_h) > 4 else 1))
    if tn_h:
        boxes.append(("Transposon", tn_h, tn_l, 1))
    if seq_h:
        boxes.append(("Sequence tracks", seq_h, seq_l, 1))

    want_legend = options.show_legend and bool(boxes)
    show_cbar = density_mappable is not None

    height = fig.get_size_inches()[1]
    fs = LEGEND_FONTSIZE
    n_cells = len(boxes) if want_legend else 0
    common = dict(frameon=True, fontsize=fs, borderpad=0.8, labelspacing=0.6,
                  columnspacing=2.2, handlelength=1.4,
                  title_fontsize=fs * TITLE_FONTSIZE_RATIO)

    # Measure the cells before reserving the band: how tall the band must be depends on
    # how many rows the cells take, and the cells must be at least as wide as the widest
    # group measures — expanding a group into a narrower cell makes matplotlib crush its
    # columns into each other.
    span = BAND_RIGHT - BAND_LEFT
    n_cols = BAND_COLS
    cell_w = (span - BAND_GAP * (n_cols - 1)) / n_cols
    cell_h_in = CBAR_ONLY_HEIGHT_IN
    if want_legend:
        n_rows_in_box = max(-(-len(handles) // ncol) for _, handles, _, ncol in boxes)
        boxes = [(title, *_pad_to_rows(handles, labels, ncol, n_rows_in_box), ncol)
                 for title, handles, labels, ncol in boxes]
        probes = [fig.legend(h, l, loc="upper left", ncol=ncol, title=t,
                             bbox_to_anchor=(0.0, 0.5), **common)
                  for t, h, l, ncol in boxes]
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        extents = [leg.get_window_extent(renderer) for leg in probes]
        w_max = max(e.width for e in extents) / fig.bbox.width
        cell_h_in = max(e.height for e in extents) / fig.dpi
        for leg in probes:
            leg.remove()
        # Two cells only if two fit; otherwise one per row, each spanning the band.
        if BAND_COLS * w_max + BAND_GAP * (BAND_COLS - 1) > span:
            n_cols = 1
        cell_w = max((span - BAND_GAP * (n_cols - 1)) / n_cols, w_max)

    # A leftover odd cell spans the row on its own, so every row shares the same left
    # and right edges.
    items: list[int | None] = list(range(n_cells)) + ([None] if show_cbar else [])
    band_rows = [items[i:i + n_cols] for i in range(0, len(items), n_cols)]
    row_w = n_cols * cell_w + BAND_GAP * (n_cols - 1)

    # Reserve a top margin for the title and a bottom band for the rows of cells.
    top_in = 0.90 if big else 0.10

    def margins(fig_h: float, n_rows: int) -> tuple[float, float]:
        """Bottom margin and band height, in inches, for ``n_rows`` rows of cells."""
        if not n_rows:
            return 0.0, 0.0
        band = n_rows * cell_h_in + (n_rows - 1) * BAND_ROW_GAP_IN + 0.10
        ruler = 0.04 * fig_h + 0.15  # pyGenomeViz hangs the ruler a height-scaled drop
        return ruler + band + 0.10, band

    bottom_in, band_in = margins(height, len(band_rows))
    if len(band_rows) > 1:
        # A second row of cells must grow the canvas, not squeeze the tracks into what
        # is left. Hold the tracks at the height they would have had with one row, then
        # settle the figure height: the ruler's drop scales with it, so iterate.
        track_in = max(0.5, height - top_in - margins(height, 1)[0])
        for _ in range(3):
            height = top_in + track_in + bottom_in
            bottom_in, band_in = margins(height, len(band_rows))
        fig.set_size_inches(fig.get_size_inches()[0], height)

    # _compress_axes scales the tracks by whatever fraction the margins leave over. On a
    # figure with few sub-tracks that fraction collapses, so grow the canvas until the
    # gene lane reaches its target. The ruler's drop scales with the height: iterate.
    lane_frac = _gene_lane_height_frac(fig)
    # `_gene_lane_height_frac` measures the whole arrow *axis*, but only `arrow_fraction`
    # of it is the drawn arrow (the rest is border pad and any histogram headroom). Pin
    # the arrow itself to MIN_GENE_LANE_RATIO so it keeps the same thickness whether or
    # not the read histogram is claiming space above it.
    lane_in = MIN_GENE_LANE_RATIO * options.track_height / max(arrow_fraction, 1e-6)
    for _ in range(4):
        needed = top_in + bottom_in + lane_in / max(lane_frac, 1e-6)
        if height >= needed:
            break
        height = needed
        bottom_in, band_in = margins(height, len(band_rows))
        fig.set_size_inches(fig.get_size_inches()[0], height)

    _compress_axes(fig, top_in, bottom_in)

    if big:
        fig.suptitle(rec.name, y=1 - 0.30 / height, fontsize=16, fontweight="bold")
        fig.text(0.5, 1 - 0.58 / height, f"{rec.accession} · {rec.length:,} bp",
                 ha="center", va="top", fontsize=10.5, color="#555555")

    if not band_rows:
        return

    band_top = _ruler_bottom_frac(fig, fallback=band_in / height) - 0.01
    cell_h = cell_h_in / height
    row_gap = BAND_ROW_GAP_IN / height

    # Place the legends first, then the colorbar cells, whose frames must be drawn to
    # match a legend's *measured* frame: matplotlib offsets a legend's frame from its
    # anchor, so a frame drawn at the nominal cell rect would sit a little high and wide.
    pending: list[tuple[float, float, float]] = []
    reference = None
    for r, row in enumerate(band_rows):
        top = band_top - r * (cell_h + row_gap)
        w = cell_w if len(row) == n_cols else row_w  # a lone cell spans the row
        x = 0.5 - row_w / 2.0  # every row centred on the same edges
        for i in row:
            if i is None:
                pending.append((x, top, w))
            else:
                title, handles, labels, ncol = boxes[i]
                leg = fig.legend(handles, labels, loc="upper left", ncol=ncol,
                                 title=title, mode="expand",
                                 bbox_to_anchor=(x, top, w, 0.0), **common)
                leg.get_frame().set_edgecolor(FRAME_EDGE_COLOR)
                reference = reference or (leg, x, top, w)
            x += w + BAND_GAP

    if pending:
        dx0 = dx1 = dy1 = 0.0
        cbar_h = cell_h
        if reference is not None:
            leg, x0, top0, w0 = reference
            fig.canvas.draw()
            e = leg.get_window_extent(fig.canvas.get_renderer())
            dx0 = e.x0 / fig.bbox.width - x0
            dx1 = e.x1 / fig.bbox.width - (x0 + w0)
            dy1 = e.y1 / fig.bbox.height - top0
            cbar_h = e.height / fig.bbox.height
        for x, top, w in pending:
            _draw_colorbar_cell(fig, density_mappable, x + dx0, top + dy1,
                                w + dx1 - dx0, cbar_h, fs)


def _ruler_bottom_frac(fig, fallback: float) -> float:
    """Figure-fraction y of the bottom of the kb scale labels (below the tracks).

    Legends are anchored just under this so they never collide with the ruler,
    which pyGenomeViz places a height-dependent distance below the lowest track.
    """
    try:
        fig.canvas.draw()
        best = None
        for ax in fig.axes:
            labels = [t for t in ax.get_xticklabels()
                      if t.get_text().strip() and "b" in t.get_text().lower()]
            if not labels:
                continue
            y0 = min(t.get_window_extent().y0 for t in labels) / fig.bbox.height
            best = y0 if best is None else min(best, y0)
        return best if best is not None else fallback
    except Exception:  # pragma: no cover - measurement is best-effort
        return fallback
