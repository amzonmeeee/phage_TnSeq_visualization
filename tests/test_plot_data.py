"""Smoke tests for sequencing-data overlays in the rendered genome map."""

from __future__ import annotations

from pathlib import Path

from phage_tnseq_viz.genome import gene_identifier, load_genome
from phage_tnseq_viz.plot import PlotOptions, render


EXAMPLES = Path(__file__).resolve().parent.parent / "examples"


def test_render_final_dataset_overlays_read_histogram_and_essentiality(tmp_path: Path):
    record = load_genome(EXAMPLES / "example_annotated.gbk")[0]
    first_gene = record.genes[0]
    output = tmp_path / "tnseq-map.svg"

    render(
        [record],
        {record.accession: []},
        PlotOptions(show_insertion_sites=False, show_insertion_density=False),
        output,
        read_counts={record.accession: {first_gene.start: 10.0}},
        gene_calls={
            (record.accession, gene_identifier(first_gene, contig=record.accession)): "Essential"
        },
    )

    svg = output.read_text(encoding="utf-8")
    assert output.exists()
    assert "#1677ff" in svg  # blue per-site read-count histogram
    assert "#800026" in svg  # essential gene-arrow fill (dark-red end of the ramp)
    assert "Tn-Seq data &amp; essentiality" in svg


def test_read_histogram_percentile_cap_marks_scale_and_renders(tmp_path: Path):
    record = load_genome(EXAMPLES / "example_annotated.gbk")[0]
    output = tmp_path / "capped.svg"

    # One hypersaturated site plus many small ones: the 90th-percentile cap must
    # sit well below the 10_000 outlier, so the scale is annotated as clipped.
    counts = {position: 5.0 for position in range(100, 300, 2)}
    counts[150] = 10_000.0

    render(
        [record],
        {record.accession: []},
        PlotOptions(
            show_insertion_sites=False,
            show_insertion_density=False,
            read_histogram_cap_percentile=90.0,
        ),
        output,
        read_counts={record.accession: counts},
    )

    svg = output.read_text(encoding="utf-8")
    assert output.exists()
    # The "≥" prefix on the read scale means the cap clipped at least one bar.
    assert "≥" in svg


def test_legacy_render_without_final_dataset_keeps_black_gene_fill(tmp_path: Path):
    record = load_genome(EXAMPLES / "example_annotated.gbk")[0]
    output = tmp_path / "previsualization.svg"

    render(
        [record],
        {record.accession: []},
        PlotOptions(show_insertion_sites=False, show_insertion_density=False),
        output,
    )

    svg = output.read_text(encoding="utf-8")
    assert "#000000" in svg
    assert "no data yet" in svg
