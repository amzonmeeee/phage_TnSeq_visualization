"""Tests for canonical final-dataset CSV input/output."""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from dataclasses import dataclass

from phage_tnseq_viz.dataset import (
    DatasetError,
    FinalSite,
    fill_missing_final_sites,
    group_counts_for_plotting,
    load_final_dataset,
    write_final_dataset,
    write_prot_table,
    write_wig,
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


@dataclass
class _Gene:
    start: int
    end: int
    strand: int
    product: str | None = None
    function: str | None = None
    locus: str | None = None


@dataclass
class _Record:
    accession: str
    genes: list


def test_write_wig_emits_transit_variablestep_with_raw_counts(tmp_path: Path):
    sites = [
        FinalSite("ctgB", 30, 0.0),
        FinalSite("ctgA", 20, 4_000_000.0),  # a hypersaturated site
        FinalSite("ctgA", 10, 7.0),
    ]

    path = write_wig(tmp_path / "out.wig", sites, comment="hello")
    text = path.read_text().splitlines()

    assert text[0] == "# hello"
    # Contigs sorted, each its own section, positions ascending within a section.
    assert text[1] == "variableStep chrom=ctgA"
    assert text[2] == "10 7"
    assert text[3] == "20 4000000"  # integer, never 4e+06
    assert text[4] == "variableStep chrom=ctgB"
    assert text[5] == "30 0"


def test_write_wig_keeps_fractional_counts_but_integers_stay_integral(tmp_path: Path):
    path = write_wig(tmp_path / "f.wig", [FinalSite("c", 1, 12.5), FinalSite("c", 2, 8.0)])
    lines = [ln for ln in path.read_text().splitlines() if not ln.startswith("variableStep")]

    assert "1 12.5" in lines
    assert "2 8" in lines  # not "8.0"


def test_wig_round_trips_through_the_tpp_parser(tmp_path: Path):
    from phage_tnseq_viz.pipeline import parse_tpp_wig

    original = [FinalSite("NC_1", 10, 5.0), FinalSite("NC_1", 25, 0.0), FinalSite("NC_1", 40, 9.0)]
    path = write_wig(tmp_path / "rt.wig", original)

    parsed = parse_tpp_wig(path)
    assert [(c.contig, c.position, c.count) for c in parsed] == [
        ("NC_1", 10, 5.0), ("NC_1", 25, 0.0), ("NC_1", 40, 9.0)
    ]


def test_write_prot_table_matches_transit_column_indices(tmp_path: Path):
    records = [
        _Record("NC_1", [
            _Gene(start=201, end=800, strand=1, product="terminase large subunit", locus="phi_001"),
            _Gene(start=996, end=1295, strand=-1, function="tail", locus="phi_002"),
            _Gene(start=1500, end=1700, strand=1),  # no product, no locus
        ])
    ]

    path = write_prot_table(tmp_path / "x.prot_table", records)
    rows = [line.split("\t") for line in path.read_text().splitlines()]

    # TRANSIT reads: description[0], start[1], end[2], strand[3], gene_name[7], orf[8].
    first = rows[0]
    assert first[0] == "terminase large subunit"
    assert first[1] == "201" and first[2] == "800" and first[3] == "+"
    assert first[7] == "-"  # no separate gene-name qualifier is captured
    assert first[8] == "phi_001"
    # aa length = (end - start + 1)//3 - 1
    assert first[4] == str((800 - 201 + 1) // 3 - 1)

    assert rows[1][3] == "-"  # reverse strand
    assert rows[1][0] == "tail"  # falls back to the PHROG function when no product
    # A gene with no locus gets a deterministic ORF id, never blank.
    assert rows[2][8] == "NC_1:1500-1700"


def test_prot_table_description_is_kept_on_one_tab_free_line(tmp_path: Path):
    records = [_Record("NC_1", [_Gene(1, 90, 1, product="messy\tname\nwith breaks", locus="g1")])]

    path = write_prot_table(tmp_path / "m.prot_table", records)
    lines = path.read_text().splitlines()

    assert len(lines) == 1  # the embedded newline did not split the row
    assert lines[0].split("\t")[0] == "messy name with breaks"
