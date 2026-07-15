"""Tests for canonical final-dataset CSV input/output."""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from phage_tnseq_viz.dataset import (
    DatasetError,
    FinalSite,
    fill_missing_final_sites,
    group_counts_for_plotting,
    load_final_dataset,
    write_final_dataset,
)
from phage_tnseq_viz.pipeline import AveragedInsertionCount, InsertionCount
from phage_tnseq_viz.dataset import final_sites_from_counts


def test_loader_accepts_aliases_coalesces_duplicates_and_default_contig(tmp_path: Path):
    path = tmp_path / "counts.csv"
    path.write_text("TA_site,count\n10,2\n10,3\n20,0\n", encoding="utf-8")

    sites = load_final_dataset(path, default_contig="phage")

    assert sites == [
        FinalSite("phage", 10, 5.0),
        FinalSite("phage", 20, 0.0),
    ]


def test_loader_rejects_missing_contig_for_multi_contig_schema(tmp_path: Path):
    path = tmp_path / "counts.csv"
    path.write_text("position,read_count\n10,2\n", encoding="utf-8")

    with pytest.raises(DatasetError, match="contig"):
        load_final_dataset(path)


def test_fill_missing_and_write_canonical_csv_with_gene_annotations(tmp_path: Path):
    sites = fill_missing_final_sites(
        [FinalSite("phage", 10, 3.0)], {"phage": [10, 20]}
    )
    path = write_final_dataset(
        tmp_path / "final_sites.csv",
        sites,
        gene_assignments={("phage", 10): (("geneA", "+"), ("geneB", "-"))},
    )

    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert [(row["position"], row["read_count"]) for row in rows] == [("10", "3"), ("20", "0")]
    assert rows[0]["gene_ids"] == "geneA;geneB"
    assert rows[0]["gene_strands"] == "+;-"


def test_pipeline_count_conversion_preserves_average_metadata_and_plot_groups():
    sites = final_sites_from_counts(
        [
            InsertionCount("phage", 10, 4),
            AveragedInsertionCount("phage", 20, 2.5, 0.5, 3),
        ]
    )

    assert sites[0].raw_read_count == 4
    assert sites[1].read_count_sd == 0.5
    assert sites[1].n_subsamples == 3
    assert group_counts_for_plotting(sites) == {"phage": {10: 4.0, 20: 2.5}}
