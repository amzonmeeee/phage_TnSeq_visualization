"""Gap statistics: is a run of consecutive empty insertion sites too long?

Saturation counts *how many* candidate sites in a gene were hit.  It says
nothing about *where* they were, so a gene whose empty sites are scattered and a
gene whose empty sites form one solid block score identically.  The second gene
is the interesting one: an uninterrupted stretch of dead sites is what a domain
that cannot tolerate insertion looks like.

This module scores that stretch.  Under a Bernoulli model where every candidate
site is hit independently with the library-wide probability ``p``, the length of
the longest run of misses in ``n`` sites follows an extreme-value (Gumbel)
distribution.  Comparing a gene's observed longest run against that distribution
gives a p-value, and a gene can be significant here while looking unremarkable
by saturation alone.

The approach is from Griffin et al. (2011), with the refinement TRANSIT later
adopted: the Gumbel location parameter is fixed by matching moments against the
*expected* longest run rather than by the cruder closed form.  That matters here,
because the expected run is computed exactly for genes with fewer than 20 sites
(Boyd's recurrence) and phage genes are mostly in that range.

Nothing in this module knows about genes or contigs; it works on plain sequences
of counts so it can be tested against published numbers directly.

References
----------
Griffin JE, Gawronski JD, DeJesus MA, Ioerger TR, Akerley BJ, Sassetti CM (2011).
High-resolution phenotypic profiling defines genes essential for mycobacterial
survival and cholesterol catabolism.  *PLoS Pathogens* 7(9):e1002251.

Schilling MF (1990).  The longest run of heads.  *College Mathematics Journal*
21(3):196-207.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
import math


EULER_GAMMA = 0.5772156649015328606

# Schilling's asymptotic expansion ends in two small oscillating correction
# terms.  Their exact n-dependence is negligible at the sizes that reach the
# asymptotic branch at all, so the reference implementation pins them to these
# constants and this port keeps the same values.
_R1 = 0.000016
_E1 = 0.01

# Above this many sites the asymptotic expansion replaces the exact recurrence,
# whose cost grows as n squared.
_EXACT_MAX_SITES = 20

# math.exp overflows past ~709; well before that the Gumbel CDF has already
# underflowed to zero, so clamping here only avoids a spurious OverflowError.
_EXP_LIMIT = 700.0


def longest_zero_run(counts: Iterable[float]) -> int:
    """Return the longest run of consecutive zero counts.

    ``counts`` must already be in genomic order: the caller decides what
    "consecutive" means, and this function does not sort.
    """

    best = 0
    current = 0
    for count in counts:
        if count > 0:
            current = 0
            continue
        current += 1
        if current > best:
            best = current
    return best


def expected_longest_run(n_sites: int, non_insertion_probability: float) -> float:
    """Expected longest run of misses over ``n_sites`` independent trials.

    Computed exactly below :data:`_EXACT_MAX_SITES` sites via the recurrence for
    ``F(n, k)``, the probability that the longest run is at most ``k``; the
    expectation then follows from ``sum(k * (F(n, k) - F(n, k - 1)))``.  Larger
    genes use Schilling's asymptotic form, which is accurate there and avoids the
    quadratic table.

    The exact branch is not an optimisation but a correctness requirement: phage
    genes commonly hold only a handful of candidate sites, where the asymptotic
    form is poor.
    """

    q = float(non_insertion_probability)
    if not 0.0 < q < 1.0:
        raise ValueError(
            f"non-insertion probability must be strictly between 0 and 1: {q!r}"
        )
    if n_sites < 0:
        raise ValueError(f"n_sites must not be negative: {n_sites!r}")
    if n_sites == 0:
        return 0.0

    n = int(n_sites)
    p = 1.0 - q

    if n < _EXACT_MAX_SITES:
        # F[m][k] = P(longest run over m trials <= k).  Entries with m <= k are
        # 1, which the initial fill already supplies.
        table = [[1.0] * (n + 1) for _ in range(n + 1)]
        for k in range(n):
            table[k + 1][k] = 1.0 - q ** (k + 1)
        for k in range(n + 1):
            weight = p * q ** (k + 1)
            for m in range(k + 2, n + 1):
                table[m][k] = table[m - 1][k] - weight * table[m - k - 2][k]
        return sum(k * (table[n][k] - table[n][k - 1]) for k in range(1, n + 1))

    log_inv_q = math.log(1.0 / q)
    return (
        math.log(n * p) / log_inv_q
        + EULER_GAMMA / log_inv_q
        - 0.5
        + _R1
        + _E1
    )


def gumbel_cdf(x: float, location: float, scale: float) -> float:
    """Gumbel CDF ``exp(-exp((location - x) / scale))``."""

    if scale <= 0.0:
        raise ValueError(f"scale must be positive: {scale!r}")
    exponent = (location - x) / scale
    if exponent > _EXP_LIMIT:
        return 0.0
    return math.exp(-math.exp(exponent))


def gap_pvalue(
    n_sites: int,
    longest_run: int,
    insertion_probability: float,
) -> float:
    """Probability of a run at least this long arising by chance.

    Returns 1.0 for the cases that carry no evidence: a gene with no candidate
    sites, a library where nothing was hit anywhere (every run is expected), and
    a fully saturated library with no run to explain.  A fully saturated library
    that nonetheless shows a run is impossible under the model and returns 0.0.
    """

    if n_sites <= 0:
        return 1.0
    if longest_run <= 0:
        return 1.0

    p = float(insertion_probability)
    if not 0.0 <= p <= 1.0:
        raise ValueError(f"insertion probability must be between 0 and 1: {p!r}")
    q = 1.0 - p
    if q >= 1.0:
        return 1.0
    if q <= 0.0:
        return 0.0

    scale = 1.0 / math.log(1.0 / q)
    location = expected_longest_run(n_sites, q) - EULER_GAMMA * scale
    return 1.0 - gumbel_cdf(longest_run, location, scale)


def benjamini_hochberg(pvalues: Sequence[float]) -> list[float]:
    """Benjamini-Hochberg adjusted p-values, in the input order.

    Values are capped at 1.0, which the raw procedure does not guarantee, and
    equal p-values receive equal q-values.
    """

    n = len(pvalues)
    if n == 0:
        return []
    for value in pvalues:
        if not 0.0 <= value <= 1.0:
            raise ValueError(f"p-values must be between 0 and 1: {value!r}")

    order = sorted(range(n), key=lambda index: pvalues[index], reverse=True)
    adjusted = [0.0] * n
    running = 1.0
    for position, index in enumerate(order):
        rank = n - position
        candidate = min(1.0, pvalues[index] * n / rank)
        running = min(running, candidate)
        adjusted[index] = running
    return adjusted


__all__ = [
    "EULER_GAMMA",
    "benjamini_hochberg",
    "expected_longest_run",
    "gap_pvalue",
    "gumbel_cdf",
    "longest_zero_run",
]
