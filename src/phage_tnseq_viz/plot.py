"""Render a linear phage genome map (pre-visualization, no sequencing data yet).

Layout per genome/contig, top to bottom within one feature track:
  * gene arrows  – black fill (essentiality placeholder), border colour = PHROG
                   functional category (phold scheme); de-novo ORFs get a grey
                   dashed border.
  * insertion track – short red ticks at every transposon insertion site.
  * (optional) GC content, GC skew, tRNA/CRISPR sub-tracks.

Built on pyGenomeViz (matplotlib) so output can be PNG or SVG, transparent, and
sized to a custom figure or a standard paper size.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

from pygenomeviz import GenomeViz

from . import colors
from .genome import GenomeRecord
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

# Sub-track heights relative to the gene-arrow lane (1.0 = same thickness as the
# arrows). The insertion ticks are matched to the arrow height on purpose so the
# arrows read as slim bars rather than a dominating block.
INSERTION_TRACK_RATIO = 1.0
DENSITY_TRACK_RATIO = 1.15
# GC content/skew are kept at their previous absolute height while the arrow lane
# was halved, so their ratio (relative to the now-thin arrow lane) is doubled.
GC_TRACK_RATIO = 2.3

# Gene-arrow lane height (pyGenomeViz feature_track_ratio). Halved from the 0.25
# default so the arrow / insertion / density lanes read as slim bars.
FEATURE_TRACK_RATIO = 0.125
# Gene arrow border line width.
GENE_EDGE_WIDTH = 4.95


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
    show_insertion_density: bool = False
    show_gc_content: bool = False
    show_gc_skew: bool = False
    show_trna: bool = False
    gc_window: int | None = None
    density_window: int | None = None
    density_cmap: str = "inferno"
    # styling
    insertion_color: str = "#e50000"
    gene_fill: str = "#000000"       # essentiality placeholder
    show_legend: bool = True
    big_title: bool = True           # genome name as a big heading on top
                                     # (False = small label beside the first row)


def _gene_edge(gene) -> tuple[str, str]:
    """Return (edgecolor, linestyle) for a gene arrow."""
    if not gene.annotated:
        return colors.UNANNOTATED_EDGE, "--"
    if gene.is_vfdb_amr:
        return colors.VFDB_AMR_COLOR, "-"
    return colors.category_color(gene.function), "-"


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
) -> Path:
    out_path = Path(out_path)

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
    straddling: list[tuple[tuple[str, int], object, str, str]] = []
    present_categories: set[str] = set()
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
                        fc=options.gene_fill,
                        ec=edge,
                        lw=GENE_EDGE_WIDTH,
                        linestyle=ls,
                    )
                else:
                    # Straddles a wrap boundary: draw the whole arrow ourselves and
                    # clip it to this row so it reads as one arrow sliced by the wrap.
                    straddling.append((key, gene, edge, ls))

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
                    f"ins:{rec.accession}:{idx}", ratio=INSERTION_TRACK_RATIO, ylim=(0, 1)
                )
            if rec.accession in density_signal:
                density_subtracks[key] = track.add_subtrack(
                    f"dens:{rec.accession}:{idx}", ratio=DENSITY_TRACK_RATIO, ylim=(0, 1)
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
        for key, gene, edge, ls in straddling:
            r0, r1 = row_window[key]
            _draw_clipped_arrow(
                row_tracks[key].ax, gene, r0, r1 - r0, head_length,
                edge, ls, options.gene_fill,
            )

    # ---- draw signals on sub-track axes (only available after plotfig) ----
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
        density_mappable=density_mappable,
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


def _default_window(length: int) -> int:
    return max(100, length // 100)


def _gc_step(window: int) -> int:
    """Sliding step for GC tracks — a small fraction of the window so the curve is
    finely sampled and smooth (phold-style) rather than coarsely polygonal."""
    return max(1, window // 10)


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


def _decorate(fig, records, options, *, present_categories, include_vfdb,
              include_denovo, transposon_label, insertion_site_count,
              density_mappable):
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
        gene_h.append(Patch(facecolor="black", edgecolor=color, linewidth=1.6))
        gene_l.append(label)
    if include_denovo:
        gene_h.append(Patch(facecolor="black", edgecolor=colors.UNANNOTATED_EDGE,
                            linewidth=1.6, linestyle="--"))
        gene_l.append("De-novo ORF (unannotated)")

    # Transposon legend (its own box).
    tn_h: list = []
    tn_l: list[str] = []
    if options.show_insertion_sites:
        motif = f" ({transposon_label})" if transposon_label else ""
        lbl = f"Insertion sites{motif} = {insertion_site_count:,}"
        tn_h.append(Line2D([0], [0], color=options.insertion_color, lw=1.6))
        tn_l.append(lbl)
    tn_h.append(Patch(facecolor="black", edgecolor="#000000", linewidth=1.0))
    tn_l.append("Gene fill = essentiality (no data yet)")

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

    # Each non-empty group becomes its own framed box. Keep the compact single-row
    # layout, but position the natural-width boxes with equal gaps so they neither
    # overlap nor look uneven.
    boxes: list[tuple[str, list, list, int]] = []
    if gene_h:
        boxes.append(("Gene function — PHROG category", gene_h, gene_l,
                      2 if len(gene_h) > 4 else 1))
    if tn_h:
        boxes.append(("Transposon", tn_h, tn_l, 1))
    if seq_h:
        boxes.append(("Sequence tracks", seq_h, seq_l, 1))

    want_legend = options.show_legend and bool(boxes)
    show_cbar = density_mappable is not None

    # Reserve a top margin for the title and a bottom margin for the legends. The
    # density colorbar goes on the right (vertical), not in the bottom band.
    height = fig.get_size_inches()[1]
    top_in = 0.90 if big else 0.10
    ruler_in = 0.04 * height + 0.15
    legend_in = 1.35 if want_legend else 0.0
    bottom_in = (ruler_in + legend_in + 0.10) if want_legend else 0.0
    _compress_axes(fig, top_in, bottom_in)

    if big:
        fig.suptitle(rec.name, y=1 - 0.30 / height, fontsize=16, fontweight="bold")
        fig.text(0.5, 1 - 0.58 / height, f"{rec.accession} · {rec.length:,} bp",
                 ha="center", va="top", fontsize=10.5, color="#555555")

    # Vertical density colorbar occupying the far right, beside the tracks.
    track_right = 1.0
    if show_cbar:
        track_right = 0.90
        for ax in fig.axes:
            p = ax.get_position()
            ax.set_position([p.x0, p.y0, max(1e-3, track_right - p.x0), p.height])
        positions = [ax.get_position() for ax in fig.axes]
        y0 = min(p.y0 for p in positions)
        y1 = max(p.y1 for p in positions)
        cax = fig.add_axes([track_right + 0.018, y0, 0.013, y1 - y0])
        cbar = fig.colorbar(density_mappable, cax=cax)
        cbar.set_label("Insertion density (sites / kb)", fontsize=8.5)
        cbar.ax.tick_params(labelsize=7.5)

    if want_legend:
        legend_top = _ruler_bottom_frac(fig, fallback=legend_in / height) - 0.01
        common = dict(frameon=True, fontsize=8, borderpad=0.8, labelspacing=0.6,
                      columnspacing=1.1, handlelength=1.4, title_fontsize=9.5)
        # Divide the bottom width into equal columns (one per box) and centre each
        # box in its column: three-way when GC tracks are shown, two-way otherwise.
        avail_left, avail_right = 0.04, track_right - 0.02
        n = len(boxes)
        for i, (title, handles, labels, ncol) in enumerate(boxes):
            cx = avail_left + (i + 0.5) * (avail_right - avail_left) / n
            leg = fig.legend(
                handles, labels, loc="upper center", bbox_to_anchor=(cx, legend_top),
                ncol=ncol, title=title, **common,
            )
            leg.get_frame().set_edgecolor("#999999")


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
