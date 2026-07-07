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

# Height of a gene arrow's shaft relative to the full track lane. A gene clipped
# at a wrap boundary (the piece without the arrowhead) is drawn as a plain
# rectangle at exactly this height so it matches the neighbouring arrow bodies
# instead of the taller full-height box pyGenomeViz would otherwise draw.
ARROW_SHAFT_RATIO = 0.5


@dataclass
class PlotOptions:
    # figure sizing
    fig_width: float = 15.0
    track_height: float = 1.6
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


def _row_label(rec: GenomeRecord, idx: int, r0: int, r1: int, n_rows: int) -> str:
    ident = f"{rec.name}\n{rec.accession} | {rec.length:,} bp"
    if n_rows == 1:
        return ident
    if idx == 0:
        return f"{ident}\n{_fmt_kb_range(r0, r1)}"
    return _fmt_kb_range(r0, r1)


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
        # Wrapped rows are stacked feature tracks; pyGenomeViz would otherwise
        # insert a full-height "link" band between them (meant for genome-to-genome
        # alignment ribbons). Shrink it to a thin inter-row gap.
        link_track_ratio=0.12,
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
            wt = compute_gc_content(rec.sequence, win, max(1, win // 2))
            gcc_signal[rec.accession] = (wt, max((abs(v) for v in wt.values), default=1.0) or 1.0)
        if options.show_gc_skew:
            win = options.gc_window or _default_window(rec.length)
            wt = compute_gc_skew(rec.sequence, win, max(1, win // 2))
            gcs_signal[rec.accession] = (wt, max((abs(v) for v in wt.values), default=1.0) or 1.0)

    # Keep references (keyed by genome + row index) so we can draw on the
    # sub-track axes after plotfig(), and remember each row's coordinate window.
    insertion_subtracks: dict[tuple[str, int], object] = {}
    density_subtracks: dict[tuple[str, int], object] = {}
    gc_content_subtracks: dict[tuple[str, int], object] = {}
    gc_skew_subtracks: dict[tuple[str, int], object] = {}
    row_window: dict[tuple[str, int], tuple[int, int]] = {}
    row_tracks: dict[tuple[str, int], object] = {}
    # Clipped gene bodies (the boundary-straddling piece without the arrowhead) to
    # draw by hand after plotfig() as shaft-height rectangles.
    truncated_bodies: list[tuple[tuple[str, int], int, int, str, str]] = []
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
                _row_label(rec, idx, r0, r1, n_rows), (r0, r1), labelsize=11
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
                # A gene straddling a wrap boundary is clipped to this row's window
                # (pyGenomeViz rejects out-of-range coordinates). Only the piece that
                # holds the gene's real head end (3' end: right for +, left for -) is
                # drawn with an arrowhead; the other piece is a flat shaft-height
                # rectangle, so the gene reads as one truncated arrow continuing
                # across rows rather than two arrows or an oversized box.
                head_in_row = g1 <= r1 if gene.strand >= 0 else g0 >= r0
                if head_in_row:
                    seg.add_feature(
                        max(g0, r0),
                        min(g1, r1),
                        gene.strand,
                        plotstyle="bigarrow",
                        arrow_shaft_ratio=ARROW_SHAFT_RATIO,
                        fc=options.gene_fill,
                        ec=edge,
                        lw=1.1,
                        linestyle=ls,
                    )
                else:
                    truncated_bodies.append((key, max(g0, r0), min(g1, r1), edge, ls))

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
                    f"ins:{rec.accession}:{idx}", ratio=0.35, ylim=(0, 1)
                )
            if rec.accession in density_signal:
                density_subtracks[key] = track.add_subtrack(
                    f"dens:{rec.accession}:{idx}", ratio=0.45, ylim=(0, 1)
                )
            if options.show_gc_content:
                gc_content_subtracks[key] = track.add_subtrack(
                    f"gcc:{rec.accession}:{idx}", ratio=0.6, ylim=(-1, 1)
                )
            if options.show_gc_skew:
                gc_skew_subtracks[key] = track.add_subtrack(
                    f"gcs:{rec.accession}:{idx}", ratio=0.6, ylim=(-1, 1)
                )

    gv.set_scale_xticks(labelsize=10)
    fig = gv.plotfig()

    # ---- truncated gene bodies: flat shaft-height rectangles at wrap cuts ----
    if truncated_bodies:
        from matplotlib.patches import Rectangle

        for key, c0, c1, edge, ls in truncated_bodies:
            track = row_tracks[key]
            ax = track.ax
            x0, x1 = track.transform_coord(c0), track.transform_coord(c1)
            y_lo, y_hi = ax.get_ylim()
            cy = (y_lo + y_hi) / 2.0
            h = ARROW_SHAFT_RATIO * (y_hi - y_lo)
            ax.add_patch(
                Rectangle(
                    (x0, cy - h / 2.0), x1 - x0, h,
                    facecolor=options.gene_fill, edgecolor=edge,
                    linewidth=1.1, linestyle=ls, joinstyle="miter", zorder=3,
                )
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

    if options.show_legend:
        _add_legend(
            fig, present_categories,
            include_vfdb=any_vfdb, include_denovo=any_denovo,
            show_insertion=options.show_insertion_sites,
            insertion_color=options.insertion_color,
            transposon_label=transposon_label,
        )

    if density_mappable is not None:
        # Reserve the colorbar strip by hand. Using fig.colorbar(ax=fig.axes)
        # steals space *unevenly* from the stacked track axes (the heatmap axis
        # ends up narrower than the gene-arrow axis), which leaves the density
        # heatmap short of the genome's 3' end. Shrinking every track axis to the
        # same right edge and placing the colorbar beside them keeps all tracks
        # left- and right-aligned.
        track_axes = list(fig.axes)
        right = max(ax.get_position().x1 for ax in track_axes)
        new_x1 = right - 0.05
        for ax in track_axes:
            bb = ax.get_position()
            ax.set_position([bb.x0, bb.y0, new_x1 - bb.x0, bb.height])
        y0 = min(ax.get_position().y0 for ax in track_axes)
        y1 = max(ax.get_position().y1 for ax in track_axes)
        cax = fig.add_axes([new_x1 + 0.012, y0, 0.012, y1 - y0])
        cbar = fig.colorbar(density_mappable, cax=cax)
        cbar.set_label("Insertion sites / kb", fontsize=8)
        cbar.ax.tick_params(labelsize=7)

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


def _add_legend(fig, categories, *, include_vfdb, include_denovo,
                show_insertion, insertion_color, transposon_label):
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    handles: list = []
    labels: list[str] = []
    for label, color in colors.legend_entries(categories, include_vfdb=include_vfdb):
        handles.append(Patch(facecolor="black", edgecolor=color, linewidth=1.6))
        labels.append(label)
    if include_denovo:
        handles.append(Patch(facecolor="black", edgecolor=colors.UNANNOTATED_EDGE,
                             linewidth=1.6, linestyle="--"))
        labels.append("De-novo ORF (unannotated)")
    if show_insertion:
        lbl = "Insertion site"
        if transposon_label:
            lbl += f" ({transposon_label})"
        handles.append(Line2D([0], [0], color=insertion_color, lw=1.5))
        labels.append(lbl)

    if not handles:
        return
    ncol = min(4, len(handles))
    fig.legend(
        handles, labels,
        loc="upper center", bbox_to_anchor=(0.5, 0.0),
        ncol=ncol, frameon=False, fontsize=9,
        title="Gene border = PHROG category   |   fill = essentiality (black = N/A)",
        title_fontsize=9,
    )
