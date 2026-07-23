"""Tests for the HMM-style confidence score.

The IQR and the per-gene confidence were cross-checked against scipy.stats.iqr
and a transcription of TRANSIT's HMM_conf ``calc_probs4``: the IQR matched to
~1e-13 and the confidence to floating-point noise over ~19,000 gene-scores.
"""

from __future__ import annotations

import math

import pytest

from phage_tnseq_viz.confidence import (
    AMBIGUOUS,
    CONFIDENT,
    LOW_CONFIDENCE_FLAG,
    interquartile_range,
    normal_pdf,
    robust_normal_params,
    score_calls,
    type7_quantile,
)


def test_type7_quantile_matches_linear_interpolation() -> None:
    values = [1.0, 2.0, 3.0, 4.0]
    assert type7_quantile(values, 0.0) == 1.0
    assert type7_quantile(values, 1.0) == 4.0
    assert type7_quantile(values, 0.5) == pytest.approx(2.5)
    assert type7_quantile(values, 0.25) == pytest.approx(1.75)
    assert type7_quantile([], 0.5) is None
    assert type7_quantile([7.0], 0.5) == 7.0


def test_interquartile_range_and_single_value() -> None:
    assert interquartile_range([1.0, 2.0, 3.0, 4.0, 5.0]) == pytest.approx(2.0)
    assert interquartile_range([5.0]) == 0.0  # undefined spread -> zero


def test_robust_params_use_median_and_floor_the_scale() -> None:
    loc, scale = robust_normal_params([10.0, 10.0, 10.0])
    assert loc == pytest.approx(10.0)
    # IQR is 0 here; the scale must not collapse, so it is floored at 1.0.
    assert scale == pytest.approx(1.0)

    # An outlier moves a mean but not the median.
    loc_robust, _ = robust_normal_params([10.0, 11.0, 12.0, 1000.0])
    assert loc_robust == pytest.approx(11.5)


def test_normal_pdf_matches_its_definition() -> None:
    assert normal_pdf(0.0, 0.0, 1.0) == pytest.approx(1.0 / math.sqrt(2 * math.pi))
    assert normal_pdf(2.0, 2.0, 3.0) == pytest.approx(1.0 / (3.0 * math.sqrt(2 * math.pi)))
    with pytest.raises(ValueError):
        normal_pdf(0.0, 0.0, 0.0)


def test_a_consistent_call_is_confident() -> None:
    """Essential genes near zero and non-essential genes high: everyone agrees."""

    items = (
        [("Essential", 0.0)] * 5
        + [("Non-essential", 200.0), ("Non-essential", 220.0), ("Non-essential", 180.0)]
    )

    scores = score_calls(items)

    assert scores is not None
    assert all(score.flag == CONFIDENT for score in scores)
    assert all(score.confidence > 0.9 for score in scores)


def test_a_call_inconsistent_with_its_count_is_flagged_low_confidence() -> None:
    """An 'Essential' gene that is actually full of reads should be doubted."""

    items = (
        [("Essential", 0.0)] * 6
        + [("Non-essential", 200.0)] * 6
        + [("Essential", 210.0)]  # called essential, but its count says otherwise
    )

    scores = score_calls(items)

    assert scores is not None
    suspect = scores[-1]
    assert suspect.label == "Essential"
    assert suspect.flag == LOW_CONFIDENCE_FLAG
    assert suspect.confidence < 0.2
    # The count clearly prefers the non-essential distribution.
    assert suspect.probabilities["Non-essential"] > suspect.probabilities["Essential"]


def test_a_borderline_count_is_flagged_ambiguous() -> None:
    """A count between two spread-out label distributions: plausible, not best.

    The groups need real spread (so their fitted scales are wide), otherwise the
    floored scale of 1.0 makes every count land unambiguously in one cluster.
    """

    essential = [0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0]  # median 15, wide IQR
    non_essential = [100.0, 120.0, 140.0, 160.0, 180.0, 200.0]  # median 150
    items = (
        [("Essential", value) for value in essential]
        + [("Non-essential", value) for value in non_essential]
        # 60 sits between the two: still plausibly essential, but non-essential
        # explains the count better.
        + [("Essential", 60.0)]
    )

    scores = score_calls(items)

    assert scores is not None
    borderline = scores[-1]
    assert borderline.flag == AMBIGUOUS
    assert borderline.confidence >= 0.2
    assert borderline.probabilities["Non-essential"] > borderline.confidence


def test_probabilities_are_normalized() -> None:
    items = [("A", 0.0), ("A", 1.0), ("B", 100.0), ("B", 110.0), ("C", 50.0), ("C", 55.0)]

    scores = score_calls(items)

    assert scores is not None
    for score in scores:
        assert sum(score.probabilities.values()) == pytest.approx(1.0)
        assert set(score.probabilities) == {"A", "B", "C"}


def test_fewer_than_two_labels_returns_none() -> None:
    assert score_calls([]) is None
    assert score_calls([("Essential", 0.0), ("Essential", 1.0)]) is None


def test_a_single_member_label_does_not_crash() -> None:
    """The scale floor keeps a one-gene label from becoming a rejecting spike."""

    scores = score_calls([("Essential", 0.0), ("Non-essential", 100.0), ("Non-essential", 120.0)])

    assert scores is not None
    assert scores[0].label == "Essential"
    assert scores[0].confidence == pytest.approx(1.0, abs=1e-6)  # 0 is unmistakably essential
    assert scores[0].flag == CONFIDENT
