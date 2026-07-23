"""Interactive, self-contained HTML genome map (optional Plotly backend).

The static :mod:`~phage_tnseq_viz.plot` renderer is ideal for print/figures, but a
long genome with thousands of insertion sites is hard to inspect at a fixed zoom.
This module renders the same information — gene arrows, per-site read counts, and
essentiality fills — as a single HTML file that pans/zooms and shows the exact
read count of each site on hover.

Plotly is an optional dependency.  It is imported lazily so the core package and
the static-image path never require it; a missing install raises a friendly
error pointing at the ``[interactive]`` extra.  The gene fill (essentiality) and
border (PHROG category) colours are shared with the static renderer so the two
outputs read as the same map.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

from . import colors
from .genome import GenomeRecord, gene_identifier
from .plot import PlotOptions, _gene_edge, _gene_fill

# Contigs are laid out left-to-right on one shared axis; this bp gap separates
# them so their arrows and a dashed boundary line never visually collide.
_CONTIG_GAP = 500

# Gene-arrow half height in the (unitless) gene row, and the head length as a
# fraction of the whole layout width, matched loosely to the static map's look.
_ARROW_HALF = 1.0
_HEAD_FRACTION = 0.012


def render_interactive_html(
    records: list[GenomeRecord],
    insertion_sites: dict[str, list[int]],
    options: PlotOptions,
    out_path: str | Path,
    transposon_label: str = "",
    *,
    read_counts: dict[str, dict[int, float]] | None = None,
    gene_calls: dict[tuple[str, str], str | None] | None = None,
) -> Path:
    """Render an interactive, offline HTML map mirroring :func:`plot.render`.

    ``read_counts`` and ``gene_calls`` use the same shapes as the static
    renderer.  The output is a single self-contained ``.html`` file (Plotly's
    library is embedded), so it opens in any browser with no network access.
    """
    try:
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
    except ImportError as exc:  # pragma: no cover - exercised only without plotly
        raise RuntimeError(
            "interactive .html output needs Plotly. Install it with "
            "`pip install \"phage-tnseq-viz[interactive]\"` or `pip install plotly`."
        ) from exc

    out_path = Path(out_path)
    read_counts = read_counts or {}
    insertion_sites = insertion_sites or {}
    show_reads = bool(options.show_read_histogram and any(read_counts.values()))
    show_ticks = bool(options.show_insertion_sites and any(insertion_sites.values()))

    total_len = sum(rec.length for rec in records) + _CONTIG_GAP * max(0, len(records) - 1)
    head = max(total_len * _HEAD_FRACTION, 1.0)
    bar_width = max(total_len / 1200.0, 3.0)

    if show_reads:
        fig = make_subplots(
            rows=2, cols=1, shared_xaxes=True, row_heights=[0.32, 0.68],
            vertical_spacing=0.06,
            subplot_titles=("Read count per insertion site", "Genes"),
        )
        gene_row = 2
    else:
        fig = make_subplots(rows=1, cols=1, subplot_titles=("Genes",))
        gene_row = 1

    seen_calls: set[str | None] = set()
    offset = 0
    boundaries: list[tuple[float, str]] = []
    for rec in records:
        _add_read_bars(
            fig, go, rec, read_counts.get(rec.accession, {}), offset, bar_width,
            options.read_histogram_color, enabled=show_reads,
        )
        _add_insertion_ticks(
            fig, go, insertion_sites.get(rec.accession, ()), offset, gene_row,
            options.insertion_color, enabled=show_ticks,
        )
        for gene in rec.genes:
            call = None
            if gene_calls is not None:
                call = gene_calls.get((rec.accession, gene_identifier(gene, contig=rec.accession)))
                seen_calls.add(call)
            _add_gene_arrow(
                fig, go, gene, rec.accession, offset, head, gene_row,
                fill=_gene_fill(rec.accession, gene, gene_calls, options.gene_fill),
                call=call,
            )
        offset += rec.length
        boundaries.append((offset, rec.accession))
        offset += _CONTIG_GAP

    _add_boundary_lines(fig, boundaries, gene_row, drawn=len(records) > 1)
    if gene_calls is not None:
        _add_essentiality_legend(fig, go, seen_calls, gene_row)

    fig.update_yaxes(
        showticklabels=False, showgrid=False, zeroline=False,
        range=[-_ARROW_HALF - 0.8, _ARROW_HALF + 0.4], fixedrange=True, row=gene_row, col=1,
    )
    if show_reads:
        fig.update_yaxes(title_text="reads", rangemode="tozero", row=1, col=1)
    fig.update_xaxes(title_text="Genomic position (bp)", row=gene_row, col=1)

    title = records[0].name if records else "Genome"
    subtitle = " · ".join(
        part for part in (records[0].accession if records else "",
                          f"{total_len:,} bp") if part
    )
    fig.update_layout(
        title=dict(text=f"{title}<br><sup>{subtitle}</sup>", x=0.5, xanchor="center"),
        bargap=0, template="plotly_white", hovermode="closest",
        legend=dict(orientation="h", yanchor="bottom", y=-0.18, xanchor="center", x=0.5,
                    title_text="Essentiality"),
        margin=dict(l=60, r=30, t=80, b=60),
    )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(str(out_path), include_plotlyjs=True, full_html=True)
    return out_path


def _add_read_bars(fig, go, rec, counts, offset, width, color, *, enabled) -> None:
    if not enabled:
        return
    points = [(position, count) for position, count in counts.items() if count > 0]
    if not points:
        return
    points.sort()
    fig.add_trace(
        go.Bar(
            x=[offset + position for position, _ in points],
            y=[count for _, count in points],
            width=width,
            marker_color=color,
            customdata=[[rec.accession, position] for position, _ in points],
            hovertemplate="%{customdata[0]}:%{customdata[1]}<br>%{y:,} reads<extra></extra>",
            showlegend=False,
        ),
        row=1, col=1,
    )


def _add_insertion_ticks(fig, go, sites, offset, gene_row, color, *, enabled) -> None:
    if not enabled or not sites:
        return
    y = -_ARROW_HALF - 0.45
    fig.add_trace(
        go.Scatter(
            x=[offset + position for position in sites],
            y=[y] * len(sites),
            mode="markers",
            marker=dict(symbol="line-ns-open", color=color, size=6, line=dict(width=0.6)),
            hovertemplate="candidate site %{x}<extra></extra>",
            showlegend=False,
        ),
        row=gene_row, col=1,
    )


def _add_gene_arrow(fig, go, gene, accession, offset, head, gene_row, *, fill, call) -> None:
    edge, linestyle = _gene_edge(gene)
    xs, ys = _arrow_xy(gene.start - 1 + offset, gene.end + offset, gene.strand, head)
    identifier = gene_identifier(gene, contig=accession)
    lines = [f"<b>{identifier}</b>", f"{accession}:{gene.start:,}–{gene.end:,} ({'+' if gene.strand >= 0 else '-'})"]
    if gene.product:
        lines.append(str(gene.product))
    lines.append(f"PHROG: {colors.PHROG_LABELS.get(colors.category_key(gene.function), 'Unknown function')}")
    if call is not None or gene.annotated:
        lines.append(f"Essentiality: {colors.essentiality_label(call)}")
    fig.add_trace(
        go.Scatter(
            x=xs, y=ys, fill="toself", mode="lines", fillcolor=fill,
            line=dict(color=edge, width=2, dash="dash" if linestyle == "--" else "solid"),
            hoveron="fills", text="<br>".join(lines), hoverinfo="text",
            showlegend=False,
        ),
        row=gene_row, col=1,
    )


def _arrow_xy(start: float, end: float, strand: int, head: float):
    """Return closed polygon coordinates for one gene arrow (home-plate shape)."""
    lo, hi, mid = -_ARROW_HALF, _ARROW_HALF, 0.0
    hl = min(head, end - start)
    if strand >= 0:  # head at the right (3') end
        xs = [start, start, end - hl, end, end - hl, start]
        ys = [lo, hi, hi, mid, lo, lo]
    else:            # head at the left (3') end
        xs = [end, end, start + hl, start, start + hl, end]
        ys = [lo, hi, hi, mid, lo, lo]
    return xs, ys


def _add_boundary_lines(fig, boundaries, gene_row, *, drawn) -> None:
    if not drawn:
        return
    # A boundary is drawn at the end of every contig except the last one.
    for end, accession in boundaries[:-1]:
        x = end + _CONTIG_GAP / 2.0
        fig.add_vline(x=x, line=dict(color="#bbbbbb", width=1, dash="dash"), row=gene_row, col=1)
    for end, accession in boundaries:
        fig.add_annotation(
            x=end - (end - 0) * 0.0, y=_ARROW_HALF + 0.35, text=accession,
            showarrow=False, xanchor="right", font=dict(size=10, color="#777777"),
            row=gene_row, col=1,
        )


def _add_essentiality_legend(fig, go, seen_calls, gene_row) -> None:
    for label, color in colors.essentiality_legend_entries(seen_calls):
        fig.add_trace(
            go.Scatter(
                x=[None], y=[None], mode="markers",
                marker=dict(size=12, color=color, line=dict(color="#555555", width=1),
                            symbol="square"),
                name=label, showlegend=True,
            ),
            row=gene_row, col=1,
        )
