"""Tests for the run-length (Griffin/Gumbel) gap statistics.

The expected-run values below were cross-checked against TRANSIT's own
implementation, which agreed to floating-point noise across 1-150 sites and
saturations from 0.05 to 0.95.
"""

from __future__ import annotations

import math

import pytest

from phage_tnseq_viz.gaps import (
    EULER_GAMMA,
    benjamini_hochberg,
    expected_longest_run,
    gap_pvalue,
    gumbel_cdf,
    longest_zero_run,
)


def test_longest_zero_run_measures_consecutive_misses() -> None:
    assert longest_zero_run([]) == 0
    assert longest_zero_run([5, 5, 5]) == 0
    assert longest_zero_run([0, 0, 0]) == 3
    assert longest_zero_run([0, 5, 0, 0, 5, 0]) == 2
    # A run that ends at the last site still counts.
    assert longest_zero_run([5, 0, 0, 0]) == 3


def test_expected_run_grows_with_sites_and_with_sparseness() -> None:
    # More trials means a longer expected run...
    assert expected_longest_run(10, 0.5) < expected_longest_run(100, 0.5)
    # ...and so does a lower chance of being hit.
    assert expected_longest_run(50, 0.2) < expected_longest_run(50, 0.8)


def test_expected_run_uses_the_exact_recurrence_for_small_genes() -> None:
    """Below 20 sites the asymptotic form is poor, and phage genes live here."""

    # Hand-checkable: with 1 site, the expected longest miss run is just q.
    assert expected_longest_run(1, 0.3) == pytest.approx(0.3)
    assert expected_longest_run(1, 0.9) == pytest.approx(0.9)

    # With 2 sites: P(run=1) = 2q(1-q), P(run=2) = q**2, so E = 2q(1-q) + 2q**2.
    q = 0.4
    assert expected_longest_run(2, q) == pytest.approx(2 * q * (1 - q) + 2 * q**2)

    # The exact and asymptotic branches meet closely at the switchover, which
    # would not hold if either branch were wrong.
    exact_side = expected_longest_run(19, 0.5)
    asymptotic_side = expected_longest_run(20, 0.5)
    assert abs(asymptotic_side - exact_side) < 0.35


def test_expected_run_rejects_degenerate_probabilities() -> None:
    for bad in (0.0, 1.0, -0.1, 1.5):
        with pytest.raises(ValueError):
            expected_longest_run(10, bad)


def test_gumbel_cdf_matches_its_definition_and_survives_underflow() -> None:
    assert gumbel_cdf(5.0, 2.0, 1.5) == pytest.approx(math.exp(-math.exp((2.0 - 5.0) / 1.5)))
    # A run far below the location would otherwise overflow math.exp.
    assert gumbel_cdf(-1e6, 0.0, 1.0) == 0.0


def test_gap_pvalue_falls_as_the_run_lengthens() -> None:
    pvalues = [gap_pvalue(40, run, 0.6) for run in range(1, 20)]

    assert pvalues == sorted(pvalues, reverse=True)
    assert pvalues[0] > 0.5  # a 1-site gap is unremarkable
    assert pvalues[-1] < 1e-6  # a 19-site gap is not


def test_gap_pvalue_accounts_for_library_saturation() -> None:
    """The same run is damning in a saturated library and ordinary in a sparse one."""

    assert gap_pvalue(30, 8, insertion_probability=0.9) < 1e-6
    assert gap_pvalue(30, 8, insertion_probability=0.2) > 0.5


def test_gap_pvalue_returns_no_evidence_for_degenerate_input() -> None:
    assert gap_pvalue(0, 0, 0.5) == 1.0  # gene with no candidate sites
    assert gap_pvalue(10, 0, 0.5) == 1.0  # no run to explain
    assert gap_pvalue(10, 5, 0.0) == 1.0  # nothing was hit anywhere
    # A run cannot occur in a fully saturated library, so observing one is
    # impossible under the model.
    assert gap_pvalue(10, 5, 1.0) == 0.0


def test_benjamini_hochberg_adjusts_and_preserves_input_order() -> None:
    pvalues = [0.001, 0.5, 0.01, 1.0]

    adjusted = benjamini_hochberg(pvalues)

    assert len(adjusted) == len(pvalues)
    assert adjusted[0] == pytest.approx(0.004)  # 0.001 * 4/1
    assert adjusted[2] == pytest.approx(0.02)  # 0.01 * 4/2
    assert all(q >= p for p, q in zip(pvalues, adjusted))
    assert all(0.0 <= q <= 1.0 for q in adjusted)


def test_benjamini_hochberg_is_monotone_and_ties_agree() -> None:
    adjusted = benjamini_hochberg([0.01, 0.01, 0.04, 0.9])

    assert adjusted[0] == adjusted[1]
    ranked = [q for _p, q in sorted(zip([0.01, 0.01, 0.04, 0.9], adjusted))]
    assert ranked == sorted(ranked)


def test_benjamini_hochberg_handles_empty_and_rejects_out_of_range() -> None:
    assert benjamini_hochberg([]) == []
    with pytest.raises(ValueError):
        benjamini_hochberg([0.5, 1.5])


def test_euler_gamma_constant() -> None:
    assert EULER_GAMMA == pytest.approx(0.5772156649, abs=1e-10)
