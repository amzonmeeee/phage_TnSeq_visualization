"""Tests for the library QC metrics.

The metric values here were cross-checked against TRANSIT's ``tnseq_stats`` and
scipy: the simple statistics agreed exactly, skewness and kurtosis to ~1e-12,
and the Pickands tail index matched TRANSIT exactly on every dataset large
enough for TRANSIT's fixed rank range.
"""

from __future__ import annotations

import random

import pytest

from phage_tnseq_viz.qc import (
    DatasetStats,
    as_row,
    central_moment_kurtosis,
    central_moment_skewness,
    dataset_stats,
    format_table,
    pickands_tail_index,
    stats_by_contig,
)


def test_basic_metrics_describe_the_count_table() -> None:
    stats = dataset_stats([0, 0, 10, 20, 30, 0], dataset="ctg")

    assert stats.dataset == "ctg"
    assert stats.sites == 6
    assert stats.hit_sites == 3
    assert stats.density == pytest.approx(0.5)
    assert stats.mean_count == pytest.approx(60 / 6)  # includes the zeros
    assert stats.nz_mean == pytest.approx(60 / 3)  # excludes them
    assert stats.nz_median == pytest.approx(20)
    assert stats.max_count == 30
    assert stats.total_counts == 60


def test_moments_are_measured_over_hit_sites_only() -> None:
    """Zeros are already reported by density; counting them again would make
    every sparse library look skewed for the wrong reason."""

    counts = [5.0, 6.0, 7.0, 100.0]
    with_zeros = dataset_stats(counts + [0.0] * 50)
    without = dataset_stats(counts)

    assert with_zeros.skewness == pytest.approx(without.skewness)
    assert with_zeros.kurtosis == pytest.approx(without.kurtosis)
    assert with_zeros.density < without.density


def test_skewness_sign_and_symmetry() -> None:
    assert central_moment_skewness([1.0, 2.0, 3.0, 4.0, 5.0]) == pytest.approx(0.0)
    assert central_moment_skewness([1.0, 1.0, 1.0, 10.0]) > 0  # long right tail
    assert central_moment_skewness([1.0, 10.0, 10.0, 10.0]) < 0


def test_moments_are_undefined_for_constant_or_tiny_samples() -> None:
    assert central_moment_skewness([7.0]) is None
    assert central_moment_kurtosis([7.0]) is None
    assert central_moment_skewness([3.0, 3.0, 3.0]) is None
    assert central_moment_kurtosis([3.0, 3.0, 3.0]) is None


def test_excess_kurtosis_is_zero_for_a_normal_sample() -> None:
    rng = random.Random(4)
    sample = [rng.gauss(0.0, 1.0) for _ in range(20000)]

    assert central_moment_kurtosis(sample) == pytest.approx(0.0, abs=0.15)


def test_pickands_index_rises_with_tail_weight() -> None:
    rng = random.Random(9)
    light = [rng.gauss(500.0, 50.0) for _ in range(2000)]
    heavy = [rng.paretovariate(0.8) * 100 for _ in range(2000)]

    assert pickands_tail_index(light) < pickands_tail_index(heavy)


def test_pickands_index_needs_enough_sites_and_says_so() -> None:
    """TRANSIT's fixed M = 10..99 scan raises IndexError below 397 sites, which
    a phage genome routinely is.  Here the scan is capped by the data instead."""

    rng = random.Random(1)
    counts = [rng.paretovariate(1.5) * 50 for _ in range(500)]

    assert pickands_tail_index(counts[:40]) is None  # too few to estimate
    assert pickands_tail_index(counts[:41]) is not None  # smallest workable set
    assert pickands_tail_index(counts[:200]) is not None  # TRANSIT would crash here


def test_pickands_index_skips_tied_order_statistics() -> None:
    """A sparse library's tail runs into the zeros, where the estimator is undefined."""

    assert pickands_tail_index([0.0] * 500) is None
    # A handful of counts over a sea of zeros still leaves no usable ranks.
    assert pickands_tail_index([900.0, 800.0, 700.0] + [0.0] * 500) is None


def test_empty_input_produces_a_zeroed_result_rather_than_an_error() -> None:
    stats = dataset_stats([], dataset="empty")

    assert stats.sites == 0
    assert stats.density == 0.0
    assert stats.skewness is None
    assert stats.pickands_tail_index is None
    assert stats.warnings  # a contig with no sites is itself worth flagging


def test_warnings_fire_on_the_published_thresholds() -> None:
    sparse = dataset_stats([0.0] * 90 + [500.0] * 10)
    assert any("saturation" in warning for warning in sparse.warnings)

    shallow = dataset_stats([1.0, 2.0, 3.0, 4.0])
    assert any("mean count at hit sites" in warning for warning in shallow.warnings)

    outlier = dataset_stats([5e6] + [100.0] * 99)
    assert any("outlier" in warning for warning in outlier.warnings)


def test_a_healthy_library_raises_no_warnings() -> None:
    rng = random.Random(2)
    counts = [0.0 if rng.random() < 0.15 else rng.uniform(80.0, 300.0) for _ in range(600)]

    assert dataset_stats(counts, dataset="healthy").warnings == ()


def test_stats_are_grouped_per_contig_in_name_order() -> None:
    rows = [("B", 5.0), ("A", 0.0), ("B", 0.0), ("A", 10.0), ("A", 20.0)]

    results = stats_by_contig(rows)

    assert [stats.dataset for stats in results] == ["A", "B"]
    assert results[0].sites == 3
    assert results[1].sites == 2
    assert results[0].total_counts == 30.0


def test_row_rendering_blanks_undefined_metrics() -> None:
    row = as_row(dataset_stats([0.0, 0.0, 5.0], dataset="ctg"))

    assert row["dataset"] == "ctg"
    assert row["density"] == "0.333"
    # One hit site cannot support a moment or a tail index.
    assert row["skewness"] == ""
    assert row["pickands_tail_index"] == ""


def test_table_renders_a_header_and_one_line_per_contig() -> None:
    lines = format_table(stats_by_contig([("A", 5.0), ("A", 0.0), ("B", 7.0)]))

    assert len(lines) == 3
    assert "contig" in lines[0] and "density" in lines[0]
    assert "A" in lines[1] and "B" in lines[2]
    assert format_table([]) == []
