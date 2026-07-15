"""End-to-end CLI tests that do not require external bioinformatics binaries."""

from __future__ import annotations

from pathlib import Path
import subprocess
import sys

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
