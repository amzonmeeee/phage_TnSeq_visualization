"""Unit tests for the non-plotting core (transposon scanning, colours, loading)."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

from phage_tnseq_viz import colors
from phage_tnseq_viz.genome import load_genome
from phage_tnseq_viz.tracks import compute_insertion_density
from phage_tnseq_viz.transposon import find_insertion_sites, resolve_transposon

EXAMPLES = Path(__file__).resolve().parent.parent / "examples"


# ---- transposon -----------------------------------------------------------

def test_mariner_preset_is_ta():
    tn = resolve_transposon("mariner")
    assert tn.motif == "TA"
    assert tn.has_preference


def test_custom_motif_uppercased():
    tn = resolve_transposon("ntan")
    assert tn.motif == "NTAN"


def test_invalid_motif_rejected():
    with pytest.raises(ValueError):
        resolve_transposon("TAXZ")


def test_random_transposon_has_no_preference():
    assert not resolve_transposon("tn5").has_preference


def test_find_ta_sites_positions():
    # TA at 0-based 2 and 6 -> 1-based 3 and 7
    seq = "GGTACCTACC"
    sites = find_insertion_sites(seq, "TA")
    assert sites == [3, 7]


def test_overlapping_motif_matches():
    # "TATA" contains TA at positions 1 and 3 (1-based)
    assert find_insertion_sites("TATA", "TA") == [1, 3]


def test_degenerate_motif_both_strands():
    # motif "GS" = G then [GC]; reverse complement scanning adds hits
    sites = find_insertion_sites("GCGC", "GC", both_strands=True)
    assert sites  # non-empty


def test_insertion_density_counts_per_kb():
    # 10 sites evenly spread in first 1000 bp of a 2000 bp genome
    sites = list(range(50, 1050, 100))  # 50,150,...,950  -> 10 sites
    track = compute_insertion_density(sites, length=2000, window=1000, step=1000)
    # first window (0-1000] holds all 10 sites -> 10/1000*1000 = 10 per kb
    assert track.values[0] == pytest.approx(10.0)
    # second window (1000-2000] holds none
    assert track.values[1] == pytest.approx(0.0)


# ---- colours --------------------------------------------------------------

def test_category_color_known():
    assert colors.category_color("tail") == "#74ee15"
    assert colors.category_color("Head and packaging") == "#ff008d"


def test_category_color_unknown_falls_back():
    assert colors.category_color(None) == colors.PHROG_COLORS["unknown function"]
    assert colors.category_color("nonsense") == colors.PHROG_COLORS["unknown function"]


# ---- genome loading -------------------------------------------------------

@pytest.fixture(scope="module", autouse=True)
def _make_examples():
    if not (EXAMPLES / "example_annotated.gbk").exists():
        subprocess.run([sys.executable, str(EXAMPLES / "make_example.py")], check=True)


def test_load_annotated_reads_functions():
    recs = load_genome(EXAMPLES / "example_annotated.gbk")
    assert len(recs) == 1
    rec = recs[0]
    assert rec.annotation_source == "genbank"
    assert rec.n_genes > 0
    assert all(g.annotated for g in rec.genes)
    assert any(g.function == "tail" for g in rec.genes)


def test_load_bare_calls_orfs():
    recs = load_genome(EXAMPLES / "example_bare.gbk")
    rec = recs[0]
    assert rec.annotation_source == "pyrodigal-gv"
    assert all(not g.annotated for g in rec.genes)


def test_bare_without_orf_finder_has_no_genes():
    recs = load_genome(EXAMPLES / "example_bare.gbk", call_orfs_if_missing=False)
    assert recs[0].n_genes == 0
