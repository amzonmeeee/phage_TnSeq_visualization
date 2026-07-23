"""Tests for TTR normalization.

The factor and its trimmed mean were cross-checked against TRANSIT's TTRNorm and
scipy.stats.trim_mean: the trimmed mean matched exactly, and the factor to
floating-point noise, over 200 random libraries.
"""

from __future__ import annotations

import random

import pytest

from phage_tnseq_viz.normalize import (
    apply_factor,
    normalize_by_contig,
    trimmed_mean,
    ttr_factor,
)


def test_trimmed_mean_matches_the_scipy_truncation_rule() -> None:
    # int(0.05 * n) is cut from each end; below 20 values nothing is trimmed.
    assert trimmed_mean(list(range(1, 20))) == pytest.approx(10.0)  # n=19, cut 0
    assert trimmed_mean(list(range(1, 21))) == pytest.approx(10.5)  # n=20, cut 1
    assert trimmed_mean(list(range(1, 41))) == pytest.approx(20.5)  # n=40, cut 2


def test_trimmed_mean_drops_the_extreme_tails() -> None:
    values = [1.0] + [10.0] * 38 + [10_000.0]  # a low and a high outlier

    # Trimming one from each end removes both, leaving a clean mean of 10.
    assert trimmed_mean(values) == pytest.approx(10.0)


def test_trimmed_mean_edge_cases() -> None:
    assert trimmed_mean([]) is None
    assert trimmed_mean([5.0]) == pytest.approx(5.0)  # nothing trimmed
    with pytest.raises(ValueError):
        trimmed_mean([1.0, 2.0], 0.5)


def test_ttr_factor_brings_the_trimmed_non_zero_mean_to_the_target() -> None:
    # density 0.5, trimmed non-zero mean 10 -> factor = 100 / (0.5 * 10) = 20.
    counts = [0.0] * 10 + [10.0] * 10
    norm = ttr_factor(counts, target=100.0)

    assert norm.density == pytest.approx(0.5)
    assert norm.trimmed_mean == pytest.approx(10.0)
    assert norm.factor == pytest.approx(20.0)


def test_ttr_scaling_makes_a_deep_and_shallow_library_comparable() -> None:
    """A library sequenced 5x deeper must land on the same scale after TTR."""

    rng = random.Random(7)
    base = [0.0 if rng.random() < 0.4 else float(rng.randint(20, 400)) for _ in range(800)]
    deep = [value * 5 for value in base]

    shallow_factor = ttr_factor(base).factor
    deep_factor = ttr_factor(deep).factor

    # The deep library gets a 5x smaller factor, cancelling its extra depth.
    assert deep_factor == pytest.approx(shallow_factor / 5, rel=1e-9)
    shallow_scaled = [apply_factor(v, shallow_factor) for v in base]
    deep_scaled = [apply_factor(v, deep_factor) for v in deep]
    assert deep_scaled == pytest.approx(shallow_scaled)


def test_ttr_factor_is_undefined_without_positive_counts() -> None:
    assert ttr_factor([0.0] * 50).factor is None
    assert ttr_factor([]).factor is None


def test_ttr_factor_rejects_a_nonpositive_target() -> None:
    with pytest.raises(ValueError):
        ttr_factor([1.0, 2.0], target=0.0)


def test_apply_factor_leaves_counts_unchanged_when_factor_is_undefined() -> None:
    assert apply_factor(37.0, None) == 37.0
    assert apply_factor(37.0, 2.0) == 74.0


def test_normalize_by_contig_scales_each_contig_independently() -> None:
    rows = (
        [("shallow", 0.0)] * 5 + [("shallow", 10.0)] * 5
        + [("deep", 0.0)] * 5 + [("deep", 1000.0)] * 5
    )

    factors = normalize_by_contig(rows, target=100.0)

    assert factors["shallow"].factor == pytest.approx(100.0 / (0.5 * 10.0))
    assert factors["deep"].factor == pytest.approx(100.0 / (0.5 * 1000.0))
    # The shallow contig is scaled up and the deep one down, toward one target.
    assert factors["shallow"].factor > 1.0 > factors["deep"].factor


def test_normalization_preserves_saturation() -> None:
    """Scaling never creates or removes a hit, so density is invariant -- which is
    why TTR cannot change a saturation-based essentiality call."""

    counts = [0.0, 0.0, 5.0, 500.0, 0.0, 12.0]
    factor = ttr_factor(counts).factor

    before = sum(1 for c in counts if c > 0)
    after = sum(1 for c in counts if apply_factor(c, factor) > 0)
    assert before == after
