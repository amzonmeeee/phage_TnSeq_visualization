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


def test_reverse_strand_insertion_point_is_motif_5prime():
    # Asymmetric motif GATA (rev-comp = TATC) in GATATATC:
    #   forward GATA at 1-based 1 -> site 1
    #   minus strand: TATC spans forward 5..8; the motif's 5' base (the G of GATA
    #   read on the - strand) is at the right end -> forward position 8.
    assert find_insertion_sites("GATATATC", "GATA") == [1, 8]
    assert find_insertion_sites("GATATATC", "GATA", both_strands=False) == [1]


def test_palindrome_single_scan_covers_both_strands():
    # TA is its own reverse complement, so one strand already finds every site.
    seq = "TAGCTAAT"
    assert (find_insertion_sites(seq, "TA", both_strands=True)
            == find_insertion_sites(seq, "TA", both_strands=False))


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


def test_infer_category_from_product():
    inf = colors.infer_category_from_product
    assert inf("terminase large subunit") == "head and packaging"
    assert inf("tail fiber protein") == "tail"
    assert inf("DNA modification methylase") == "dna, rna and nucleotide metabolism"
    assert inf("phage integrase") == "integration and excision"
    assert inf("endolysin") == "lysis"
    assert inf("hypothetical protein") == "unknown function"
    assert inf("head-tail adaptor") == "connector"
    assert inf(None) is None
    assert inf("some totally novel widget") is None


# ---- line wrapping --------------------------------------------------------

def test_row_windows_wraps_and_balances():
    from phage_tnseq_viz.plot import _row_windows

    # default 20 kb wrap on a 39,422 bp phage -> 2 evenly balanced rows
    rows = _row_windows(39422, 20.0, None)
    assert len(rows) == 2
    assert rows[0][0] == 0 and rows[-1][1] == 39422
    # contiguous coverage, no gaps or overlaps
    for a, b in zip(rows, rows[1:]):
        assert a[1] == b[0]
    # balanced: row widths differ by at most 1 bp
    widths = [e - s for s, e in rows]
    assert max(widths) - min(widths) <= 1


def test_row_windows_single_line_cases():
    from phage_tnseq_viz.plot import _row_windows

    # short genome stays on one line under the default wrap width
    assert _row_windows(5000, 20.0, None) == [(0, 5000)]
    # wrap disabled -> single line regardless of length
    assert _row_windows(60000, 0.0, None) == [(0, 60000)]


def test_row_windows_force_rows_overrides():
    from phage_tnseq_viz.plot import _row_windows

    rows = _row_windows(30000, 20.0, 3)
    assert len(rows) == 3
    assert rows[0][0] == 0 and rows[-1][1] == 30000


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
