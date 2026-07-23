"""TTR normalization for per-site Tn-Seq counts.

Sequencing depth is arbitrary: the same library sequenced twice as deep has
twice the counts everywhere, and two phage libraries prepared side by side rarely
match. Comparing their raw counts, or drawing them on the same scale, is
therefore meaningless until they are put on a common footing.

The package already offers one way to do that -- SeqKit subsampling every library
down to a shared read depth -- but subsampling discards data and injects sampling
noise. TTR (Trimmed Total Reads) instead rescales: it finds the single factor
that brings a library's typical non-zero count to a fixed target, and multiplies
every count by it. Nothing is thrown away, and a library run on its own reaches
the same target as any other, so the results are comparable without the two ever
being processed together.

The factor is ``target / (density * trimmed_mean_of_hit_sites)``, where the
trimmed mean drops the top and bottom 5% of non-zero counts before averaging.
That trimming is the whole point of the "T": a Tn-Seq count distribution has a
heavy upper tail (a few hypersaturated sites), and an ordinary mean would let
those few sites set the scale for the whole library.

This is TRANSIT's default normalization, and the formula here matches its
implementation. Two properties are worth stating plainly, because they bound
what TTR can and cannot do:

* It rescales a whole library by one number, so it never changes any *within*
  -library comparison: saturation is untouched (zero stays zero), and every
  essentiality call is unchanged because all counts move together. Its effect is
  entirely on *cross*-library comparability and on the map's read-count scale.
* It assumes depth is the only thing separating libraries. It does not correct a
  position-dependent coverage gradient or a genuinely different insertion-site
  distribution.

Everything here is pure Python, like the other analysis modules, so it adds no
dependency.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass
import math


TTR_TARGET = 100.0
TTR_TRIM = 0.05


@dataclass(frozen=True)
class ContigNormalization:
    """The TTR factor applied to one contig, with the inputs that produced it.

    ``factor`` is ``None`` when no factor could be computed -- a contig with no
    hit sites, or one whose counts all fall inside the trimmed region -- in which
    case the counts are left unchanged.
    """

    contig: str
    factor: float | None
    density: float
    trimmed_mean: float | None
    target: float


def trimmed_mean(values: Sequence[float], proportiontocut: float = TTR_TRIM) -> float | None:
    """Mean after removing a proportion of the largest and smallest values.

    This reproduces ``scipy.stats.trim_mean``: the cut count is
    ``int(proportiontocut * n)`` at each end, using truncation, so a small sample
    may have nothing trimmed. Returns ``None`` when the trimmed region would be
    empty (which for the 5% default needs at least one value).
    """

    if not 0.0 <= proportiontocut < 0.5:
        raise ValueError(f"proportiontocut must be in [0, 0.5): {proportiontocut!r}")
    n = len(values)
    if n == 0:
        return None
    lowercut = int(proportiontocut * n)
    uppercut = n - lowercut
    if uppercut <= lowercut:
        return None
    ordered = sorted(float(value) for value in values)
    kept = ordered[lowercut:uppercut]
    return sum(kept) / len(kept)


def ttr_factor(
    counts: Sequence[float],
    *,
    target: float = TTR_TARGET,
    trim: float = TTR_TRIM,
) -> ContigNormalization:
    """Compute the TTR scaling factor for one contig's complete count table.

    ``counts`` must include zero-count sites: density is measured over every
    candidate site, so leaving the zeros out would inflate it and shrink the
    factor. The contig label is filled in by the caller via
    :func:`normalize_by_contig`; here it is left blank.
    """

    if target <= 0.0:
        raise ValueError(f"target must be positive: {target!r}")
    values = [float(count) for count in counts]
    n = len(values)
    density = (sum(1 for value in values if value > 0.0) / n) if n else 0.0
    positive = [value for value in values if value > 0.0]
    mu = trimmed_mean(positive, trim)

    if mu is None or mu <= 0.0 or density <= 0.0:
        factor: float | None = None
    else:
        factor = target / (density * mu)
    return ContigNormalization(
        contig="",
        factor=factor,
        density=density,
        trimmed_mean=mu,
        target=target,
    )


def apply_factor(count: float, factor: float | None) -> float:
    """Scale one count, leaving it unchanged when the factor is undefined."""

    if factor is None:
        return float(count)
    return float(count) * factor


def normalize_by_contig(
    rows: Iterable[tuple[str, float]],
    *,
    target: float = TTR_TARGET,
    trim: float = TTR_TRIM,
) -> dict[str, ContigNormalization]:
    """Compute a separate TTR factor for each contig's counts.

    Factors are per contig for the same reason read thresholds are elsewhere in
    this package: phage contigs are not sequenced to comparable depths, so
    pooling them would let a deep contig set the scale for a shallow one. Each
    contig is put on the common target independently.
    """

    grouped: dict[str, list[float]] = {}
    for contig, count in rows:
        grouped.setdefault(contig, []).append(float(count))
    result: dict[str, ContigNormalization] = {}
    for contig, counts in grouped.items():
        normalization = ttr_factor(counts, target=target, trim=trim)
        result[contig] = ContigNormalization(
            contig=contig,
            factor=normalization.factor,
            density=normalization.density,
            trimmed_mean=normalization.trimmed_mean,
            target=normalization.target,
        )
    return result


__all__ = [
    "TTR_TARGET",
    "TTR_TRIM",
    "ContigNormalization",
    "apply_factor",
    "normalize_by_contig",
    "trimmed_mean",
    "ttr_factor",
]
