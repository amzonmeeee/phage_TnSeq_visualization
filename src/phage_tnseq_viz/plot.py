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

from dataclasses import dataclass
from pathlib import Path

from pygenomeviz import GenomeViz

from . import colors
from .genome import GenomeRecord
from .tracks import compute_gc_content, compute_gc_skew, compute_insertion_density

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


@dataclass
class PlotOptions:
    # figure sizing
    fig_width: float = 15.0
    track_height: float = 1.6
    paper: str | None = None
    portrait: bool = False
    fit_page: bool = False
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
    )

    # Keep references so we can draw on sub-track axes after plotfig().
    insertion_subtracks: dict[str, object] = {}
    density_subtracks: dict[str, object] = {}
    gc_content_subtracks: dict[str, object] = {}
    gc_skew_subtracks: dict[str, object] = {}
    present_categories: set[str] = set()
    any_denovo = False
    any_vfdb = False

    for rec in records:
        label = f"{rec.name}\n{rec.accession} | {rec.length:,} bp"
        track = gv.add_feature_track(label, rec.length, labelsize=11)
        seg = track.get_segment()

        for gene in rec.genes:
            edge, ls = _gene_edge(gene)
            if gene.annotated and not gene.is_vfdb_amr:
                present_categories.add(colors.category_key(gene.function))
            any_denovo = any_denovo or (not gene.annotated)
            any_vfdb = any_vfdb or gene.is_vfdb_amr
            seg.add_feature(
                gene.start - 1,
                gene.end,
                gene.strand,
                plotstyle="arrow",
                fc=options.gene_fill,
                ec=edge,
                lw=1.1,
                linestyle=ls,
            )

        # optional tRNA / CRISPR as thin boxes on the main lane
        if options.show_trna:
            for f in rec.noncds:
                seg.add_feature(
                    f.start - 1, f.end, f.strand,
                    plotstyle="rbox", fc=TRNA_COLOR, ec=TRNA_COLOR, lw=0.5,
                )

        if options.show_insertion_sites and insertion_sites.get(rec.accession):
            insertion_subtracks[rec.accession] = track.add_subtrack(
                f"ins:{rec.accession}", ratio=0.35, ylim=(0, 1)
            )
        if options.show_insertion_density and insertion_sites.get(rec.accession):
            density_subtracks[rec.accession] = track.add_subtrack(
                f"dens:{rec.accession}", ratio=0.45, ylim=(0, 1)
            )
        if options.show_gc_content:
            gc_content_subtracks[rec.accession] = track.add_subtrack(
                f"gcc:{rec.accession}", ratio=0.6, ylim=(-1, 1)
            )
        if options.show_gc_skew:
            gc_skew_subtracks[rec.accession] = track.add_subtrack(
                f"gcs:{rec.accession}", ratio=0.6, ylim=(-1, 1)
            )

    gv.set_scale_xticks(labelsize=10)
    fig = gv.plotfig()

    # ---- draw signals on sub-track axes (only available after plotfig) ----
    density_mappable = None
    for rec in records:
        sub = insertion_subtracks.get(rec.accession)
        if sub is not None:
            xs = [sub.transform_coord(p) for p in insertion_sites.get(rec.accession, [])]
            if xs:
                sub.ax.vlines(xs, 0.0, 1.0, color=options.insertion_color, lw=0.4)

        sub = density_subtracks.get(rec.accession)
        if sub is not None:
            m = _draw_density(
                sub, rec, insertion_sites.get(rec.accession, []),
                window=options.density_window, cmap=options.density_cmap,
            )
            density_mappable = density_mappable or m

        sub = gc_content_subtracks.get(rec.accession)
        if sub is not None:
            _draw_window(sub, rec, kind="content", window=options.gc_window)

        sub = gc_skew_subtracks.get(rec.accession)
        if sub is not None:
            _draw_window(sub, rec, kind="skew", window=options.gc_window)

    if options.show_legend:
        _add_legend(
            fig, present_categories,
            include_vfdb=any_vfdb, include_denovo=any_denovo,
            show_insertion=options.show_insertion_sites,
            insertion_color=options.insertion_color,
            transposon_label=transposon_label,
        )

    if density_mappable is not None:
        cbar = fig.colorbar(
            density_mappable, ax=fig.axes, location="right",
            fraction=0.012, pad=0.006, aspect=30,
        )
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


def _draw_window(sub, rec: GenomeRecord, *, kind: str, window: int | None) -> None:
    win = window or _default_window(rec.length)
    step = max(1, win // 2)
    if kind == "content":
        wt = compute_gc_content(rec.sequence, win, step)
        pos_c, neg_c = GC_POSITIVE, GC_NEGATIVE
    else:
        wt = compute_gc_skew(rec.sequence, win, step)
        pos_c, neg_c = SKEW_POSITIVE, SKEW_NEGATIVE
    if not wt.values:
        return
    clamped = [min(max(p, 0), rec.length) for p in wt.positions]
    xs = [sub.transform_coord(p) for p in clamped]
    ys = wt.values
    vmax = max(abs(v) for v in ys) or 1.0
    ys = [v / vmax for v in ys]  # scale into (-1, 1)
    sub.ax.fill_between(xs, ys, 0, where=[v >= 0 for v in ys], color=pos_c, lw=0)
    sub.ax.fill_between(xs, ys, 0, where=[v < 0 for v in ys], color=neg_c, lw=0)
    sub.ax.axhline(0, color="#00000055", lw=0.4)


def _draw_density(sub, rec: GenomeRecord, sites, *, window, cmap):
    """Draw an insertion-density heat strip; return the mappable for a colorbar."""
    import numpy as np

    win = window or max(200, rec.length // 100)
    step = max(1, win // 2)
    wt = compute_insertion_density(sites, rec.length, win, step)
    if not wt.values:
        return None
    data = np.array(wt.values).reshape(1, -1)
    x0 = sub.transform_coord(0)
    x1 = sub.transform_coord(rec.length)
    im = sub.ax.imshow(
        data, extent=(x0, x1, 0, 1), aspect="auto",
        cmap=cmap, origin="lower", interpolation="nearest",
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
