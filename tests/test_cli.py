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
    assert main(common + ["-o", str(raw_dir / "m.svg"), "--csv-dir", str(raw_dir)]) == 0
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
