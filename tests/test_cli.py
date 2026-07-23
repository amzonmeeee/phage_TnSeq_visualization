"""End-to-end CLI tests that do not require external bioinformatics binaries."""

from __future__ import annotations

from pathlib import Path
import subprocess
import sys

import pytest

from phage_tnseq_viz.cli import main


EXAMPLES = Path(__file__).resolve().parent.parent / "examples"


def test_example_generator_creates_a_complete_final_dataset(tmp_path: Path):
    subprocess.run(
        [sys.executable, str(EXAMPLES / "make_example.py"), "--output-dir", str(tmp_path)],
        check=True,
    )

    assert (tmp_path / "example_annotated.gbk").exists()
    final_csv = tmp_path / "example_final_sites.csv"
    assert final_csv.exists()
    assert "read_count" in final_csv.read_text(encoding="utf-8").splitlines()[0]


def test_plot_final_dataset_bypasses_processing_and_writes_map_and_csvs(tmp_path: Path):
    counts = tmp_path / "final_counts.csv"
    # Single-contig shorthand: the CLI supplies the GenBank accession as default.
    counts.write_text("TA_site,read_count\n300,8\n500,0\n700,4\n", encoding="utf-8")
    image = tmp_path / "map.svg"

    rc = main(
        [
            "plot", str(EXAMPLES / "example_annotated.gbk"),
            "--final-dataset", str(counts),
            "--no-insertion-sites", "--no-insertion-density",
            "--output", str(image), "--csv-dir", str(tmp_path),
        ]
    )

    assert rc == 0
    assert image.exists()
    assert (tmp_path / "map_sites.csv").exists()
    assert (tmp_path / "map_gene_essentiality.csv").exists()


def test_html_extension_writes_a_self_contained_interactive_map(tmp_path: Path):
    pytest.importorskip("plotly")
    counts = tmp_path / "final_counts.csv"
    counts.write_text("TA_site,read_count\n300,8\n500,0\n700,4\n", encoding="utf-8")
    out = tmp_path / "map.html"

    rc = main(
        [
            "plot", str(EXAMPLES / "example_annotated.gbk"),
            "--final-dataset", str(counts),
            "--no-insertion-sites", "--no-insertion-density",
            "-o", str(out), "--csv-dir", str(tmp_path),
        ]
    )

    assert rc == 0
    assert out.exists()
    html = out.read_text(encoding="utf-8")
    # A real, offline Plotly document: the draw call plus an inlined library.
    assert "Plotly.newPlot" in html
    assert "reads" in html  # per-site read-count hover data is present


def test_quiet_suppresses_progress_but_verbose_and_errors_still_speak(tmp_path: Path, capsys):
    counts = tmp_path / "final_counts.csv"
    counts.write_text("TA_site,read_count\n300,8\n700,4\n", encoding="utf-8")
    common = [
        "plot", str(EXAMPLES / "example_annotated.gbk"),
        "--final-dataset", str(counts),
        "--no-insertion-sites", "--no-insertion-density",
        "--csv-dir", str(tmp_path),
    ]

    assert main(common + ["-o", str(tmp_path / "normal.png")]) == 0
    assert "Done. Wrote" in capsys.readouterr().out

    assert main(common + ["-o", str(tmp_path / "quiet.png"), "--quiet"]) == 0
    assert capsys.readouterr().out.strip() == ""

    # An error is emitted regardless of --quiet, on stderr.
    assert main(common + ["-o", str(tmp_path / "bad.png"), "--quiet", "--read-histogram-cap", "999"]) == 1
    captured = capsys.readouterr()
    assert captured.out.strip() == ""
    assert "error:" in captured.err


def test_process_can_skip_every_external_stage_and_still_records_manifest(tmp_path: Path):
    reference = EXAMPLES / "example_annotated.gbk"
    rc = main(
        [
            "process", str(reference), "--reads1", str(reference),
            "--skip-fastqc", "--skip-fastp", "--skip-tpp",
            "--output-dir", str(tmp_path / "run"),
        ]
    )

    assert rc == 0
    assert (tmp_path / "run" / "reference.fasta").exists()
    assert (tmp_path / "run" / "processing_manifest.json").exists()


def _final_calls(path: Path) -> dict[str, str]:
    import csv
    with path.open() as handle:
        return {row["gene_id"]: row["final_call"] for row in csv.DictReader(handle)}


def test_ttr_normalization_preserves_calls_and_reports_raw_qc(tmp_path: Path):
    """The plot path with --normalize ttr must: rescale read_count, keep the raw
    count in raw_read_count, leave every essentiality call unchanged, and still
    compute QC on the raw counts."""

    import csv

    reference = EXAMPLES / "example_annotated.gbk"
    dataset = EXAMPLES / "example_final_sites.csv"
    if not dataset.exists():  # generated on demand, mirroring the demo
        subprocess.run(
            [sys.executable, str(EXAMPLES / "make_example.py"), "--output-dir", str(tmp_path / "gen")],
            check=True,
        )
        reference = tmp_path / "gen" / "example_annotated.gbk"
        dataset = tmp_path / "gen" / "example_final_sites.csv"

    common = [
        "plot", str(reference), "--final-dataset", str(dataset),
        "--no-insertion-sites", "--no-insertion-density",
    ]
    raw_dir = tmp_path / "raw"
    ttr_dir = tmp_path / "ttr"
    # Pin the baseline to none; the default is now auto (which normalizes on the
    # plot path), so an un-normalized reference must be requested explicitly.
    assert main(common + ["-o", str(raw_dir / "m.svg"), "--csv-dir", str(raw_dir),
                          "--normalize", "none"]) == 0
    assert main(common + ["-o", str(ttr_dir / "m.svg"), "--csv-dir", str(ttr_dir),
                          "--normalize", "ttr"]) == 0

    # 1. Within-library invariance: calls are identical.
    assert _final_calls(raw_dir / "m_gene_essentiality.csv") == _final_calls(
        ttr_dir / "m_gene_essentiality.csv"
    )

    # 2. read_count is rescaled, raw_read_count keeps the original.
    with (ttr_dir / "m_sites.csv").open() as handle:
        rows = [row for row in csv.DictReader(handle) if float(row["read_count"]) > 0]
    assert rows
    assert all(row["raw_read_count"] for row in rows)
    assert any(
        abs(float(row["read_count"]) - float(row["raw_read_count"])) > 1e-6 for row in rows
    )

    # 3. QC is computed on raw counts, so its max matches the raw table, not the
    #    normalized one.
    with (ttr_dir / "m_qc.csv").open() as handle:
        qc_row = next(csv.DictReader(handle))
    with (raw_dir / "m_sites.csv").open() as handle:
        raw_max = max(float(row["read_count"]) for row in csv.DictReader(handle))
    assert float(qc_row["max_count"]) == pytest.approx(raw_max)


def test_plot_emits_wig_and_prot_table_and_accepts_wig_input(tmp_path: Path):
    """The TRANSIT hand-off: raw WIG + prot_table out, and WIG back in."""

    import csv

    counts = tmp_path / "final.csv"
    counts.write_text("TA_site,read_count\n300,8\n500,0\n700,4\n", encoding="utf-8")
    out_dir = tmp_path / "out"

    rc = main([
        "plot", str(EXAMPLES / "example_annotated.gbk"),
        "--final-dataset", str(counts),
        "--no-insertion-sites", "--no-insertion-density",
        "-o", str(out_dir / "m.svg"), "--csv-dir", str(out_dir),
    ])
    assert rc == 0

    wig = out_dir / "m.wig"
    prot_table = out_dir / "m.prot_table"
    assert wig.exists() and prot_table.exists()
    assert "variableStep chrom=" in wig.read_text()
    # prot_table rows have the ORF at index 8, the way TRANSIT reads it.
    first_gene = prot_table.read_text().splitlines()[0].split("\t")
    assert len(first_gene) >= 9 and first_gene[8]

    # Feed the emitted WIG straight back as the dataset; it must be accepted and
    # reproduce the same gene calls as the CSV it came from.
    from_csv = out_dir / "m_gene_essentiality.csv"
    rt_dir = tmp_path / "rt"
    rc = main([
        "plot", str(EXAMPLES / "example_annotated.gbk"),
        "--final-dataset", str(wig),
        "--no-insertion-sites", "--no-insertion-density",
        "-o", str(rt_dir / "m.svg"), "--csv-dir", str(rt_dir),
    ])
    assert rc == 0

    def calls(path: Path) -> dict[str, str]:
        with path.open() as handle:
            return {row["gene_id"]: row["final_call"] for row in csv.DictReader(handle)}

    assert calls(from_csv) == calls(rt_dir / "m_gene_essentiality.csv")


def test_ttr_normalization_leaves_the_transit_wig_raw(tmp_path: Path):
    """The WIG must carry raw counts even under --normalize ttr, because TRANSIT
    normalizes itself and would otherwise double-normalize."""

    counts = tmp_path / "final.csv"
    counts.write_text("TA_site,read_count\n300,8\n500,0\n700,40\n900,4\n", encoding="utf-8")
    out_dir = tmp_path / "out"

    rc = main([
        "plot", str(EXAMPLES / "example_annotated.gbk"),
        "--final-dataset", str(counts),
        "--no-insertion-sites", "--no-insertion-density",
        "-o", str(out_dir / "m.svg"), "--csv-dir", str(out_dir),
        "--normalize", "ttr",
    ])
    assert rc == 0

    # The WIG holds the raw inputs; TTR scaled the CSV's read_count to something else.
    wig_counts = {}
    for line in (out_dir / "m.wig").read_text().splitlines():
        if line.startswith(("#", "variableStep")):
            continue
        pos, value = line.split()
        wig_counts[int(pos)] = float(value)
    assert wig_counts[300] == 8.0
    assert wig_counts[700] == 40.0

    import csv
    with (out_dir / "m_sites.csv").open() as handle:
        normalized = {int(r["position"]): float(r["read_count"]) for r in csv.DictReader(handle)}
    # Normalization actually changed the CSV values, so the WIG staying raw is meaningful.
    assert normalized[700] != 40.0


def test_normalize_auto_resolves_by_whether_counts_came_from_tpp():
    from phage_tnseq_viz.cli import _resolve_normalize

    # The plot path (own data, not via TPP): auto normalizes here.
    assert _resolve_normalize("auto", from_tpp=False) == "ttr"
    # The process path (counts from TPP -> TRANSIT): auto defers.
    assert _resolve_normalize("auto", from_tpp=True) == "none"
    # An explicit choice always wins over auto, on either path.
    assert _resolve_normalize("none", from_tpp=False) == "none"
    assert _resolve_normalize("ttr", from_tpp=True) == "ttr"


def test_normalize_none_keeps_the_plot_counts_raw(tmp_path: Path):
    """auto normalizes the plot path, so --normalize none is how you opt out."""

    import csv

    counts = tmp_path / "final.csv"
    counts.write_text("TA_site,read_count\n300,8\n500,0\n700,40\n", encoding="utf-8")
    out = tmp_path / "out"

    assert main([
        "plot", str(EXAMPLES / "example_annotated.gbk"), "--final-dataset", str(counts),
        "--no-insertion-sites", "--no-insertion-density",
        "-o", str(out / "m.svg"), "--csv-dir", str(out), "--normalize", "none",
    ]) == 0

    with (out / "m_sites.csv").open() as handle:
        by_pos = {int(r["position"]): r for r in csv.DictReader(handle)}
    # Untouched: read_count equals the input and no raw_read_count shadow is set.
    assert float(by_pos[700]["read_count"]) == 40.0
    assert by_pos[700]["raw_read_count"] == ""
