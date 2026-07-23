"""Confidence scoring for essentiality calls.

An essentiality call is made from where the insertions fall, but a gene's overall
insertion *level* is a second, independent check on that call: essential genes
have counts near zero, non-essential genes have high counts, and the fitness
categories sit in between. When a gene's call disagrees with its own count -- a
gene called essential that is full of reads, or called non-essential with none --
the call deserves a second look.

This module quantifies that agreement without re-deriving the calls. For each
label actually assigned, it fits a robust Normal distribution to the mean counts
of the genes carrying that label, then asks, for every gene, how well its own
count fits its label's distribution relative to the others. The result is a
confidence in [0, 1] and a flag when the count points elsewhere.

The method is TRANSIT's HMM confidence score (its ``HMM_conf.py``), lifted off
the HMM: it is a post-processing step over a set of labelled genes and their mean
counts, so it applies just as well to this package's rule-based calls or a custom
classifier's labels. The robust choices that make it survive small groups are
kept intact -- a median for location and an IQR-based scale floored at 1.0 -- and
matter far more here than in TRANSIT, because a phage genome offers tens of genes
per label, not thousands.

Everything here is pure Python, like the other analysis modules.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
import math


# Below this normalized probability a call is inconsistent with its own count.
LOW_CONFIDENCE = 0.2

# Converts an interquartile range to a Normal-consistent standard deviation
# (1 / 1.349). TRANSIT's source carries a transposed 0.7314; the floor below
# dominates for the small groups where the difference could matter.
_IQR_TO_SIGMA = 0.7413

# The scale is never allowed below this, so a label whose genes share almost the
# same count cannot collapse to a spike that rejects every other gene.
_SIGMA_FLOOR = 1.0

CONFIDENT = ""
AMBIGUOUS = "ambiguous"
LOW_CONFIDENCE_FLAG = "low-confidence"


@dataclass(frozen=True)
class ConfidenceScore:
    """How well one gene's count agrees with its essentiality label.

    ``confidence`` is the normalized probability of the assigned label given the
    gene's mean count. ``flag`` is empty when the label is the most probable
    explanation of the count, ``"low-confidence"`` when the probability is below
    :data:`LOW_CONFIDENCE` (ignore the call), and ``"ambiguous"`` when the label
    is plausible but another label fits the count better.

    ``probabilities`` is the full normalized distribution over labels, so a
    reader can see which call the count would have preferred.
    """

    label: str
    confidence: float
    flag: str
    probabilities: dict[str, float]


def type7_quantile(values: Sequence[float], probability: float) -> float | None:
    """R/NumPy type-7 quantile (linear interpolation), matching ``scipy.stats``."""

    if not values:
        return None
    ordered = sorted(float(value) for value in values)
    if len(ordered) == 1:
        return ordered[0]
    h = (len(ordered) - 1) * probability
    low = math.floor(h)
    high = math.ceil(h)
    if low == high:
        return ordered[low]
    return ordered[low] + (h - low) * (ordered[high] - ordered[low])


def interquartile_range(values: Sequence[float]) -> float:
    """Q3 - Q1 using the same interpolation scheme as ``scipy.stats.iqr``."""

    if len(values) < 2:
        return 0.0
    q1 = type7_quantile(values, 0.25)
    q3 = type7_quantile(values, 0.75)
    return float(q3) - float(q1)


def robust_normal_params(values: Sequence[float]) -> tuple[float, float]:
    """Return (location, scale) as a median and an IQR-derived, floored sigma."""

    if not values:
        raise ValueError("cannot fit a distribution to no values")
    location = type7_quantile(values, 0.50)
    scale = max(_SIGMA_FLOOR, _IQR_TO_SIGMA * interquartile_range(values))
    return float(location), scale


def normal_pdf(x: float, location: float, scale: float) -> float:
    """Normal probability density; ``scale`` must be positive."""

    return math.exp(normal_logpdf(x, location, scale))


def normal_logpdf(x: float, location: float, scale: float) -> float:
    """Log of the Normal density; ``scale`` must be positive.

    Scoring works in log space because the scale is floored at 1.0, so a count
    far from every label's cluster gives densities that all underflow to zero at
    once -- which happens with the tight count clusters of a small phage genome,
    though not with the thousands of genes TRANSIT was built for. Normalizing the
    logs recovers the same probabilities without the zero-division.
    """

    if scale <= 0.0:
        raise ValueError(f"scale must be positive: {scale!r}")
    z = (x - location) / scale
    return -0.5 * z * z - math.log(scale) - 0.5 * math.log(2.0 * math.pi)


def score_calls(
    items: Sequence[tuple[str, float]],
) -> list[ConfidenceScore] | None:
    """Score a set of labelled genes by how well each fits its label's counts.

    ``items`` pairs each gene's assigned label with its mean count over all
    candidate sites (zeros included). The returned list aligns with ``items``.

    Returns ``None`` when the exercise is not meaningful: fewer than two distinct
    labels leaves nothing to compare against, so every gene would be trivially
    certain of the only option.
    """

    labels_present = {label for label, _mean in items}
    if len(labels_present) < 2:
        return None

    params = {
        label: robust_normal_params(
            [mean for other, mean in items if other == label]
        )
        for label in labels_present
    }

    scores: list[ConfidenceScore] = []
    for label, mean in items:
        log_densities = {
            candidate: normal_logpdf(mean, location, scale)
            for candidate, (location, scale) in params.items()
        }
        # Stable softmax: shifting the logs by their max cancels in the ratio and
        # guarantees at least one non-underflowing term, so the total is positive.
        offset = max(log_densities.values())
        weights = {name: math.exp(value - offset) for name, value in log_densities.items()}
        total = sum(weights.values())
        probabilities = {name: weight / total for name, weight in weights.items()}
        confidence = probabilities[label]
        best = max(probabilities.values())

        if confidence < LOW_CONFIDENCE:
            flag = LOW_CONFIDENCE_FLAG
        elif confidence != best:
            flag = AMBIGUOUS
        else:
            flag = CONFIDENT
        scores.append(
            ConfidenceScore(
                label=label,
                confidence=confidence,
                flag=flag,
                probabilities=probabilities,
            )
        )
    return scores


__all__ = [
    "AMBIGUOUS",
    "CONFIDENT",
    "LOW_CONFIDENCE",
    "LOW_CONFIDENCE_FLAG",
    "ConfidenceScore",
    "interquartile_range",
    "normal_logpdf",
    "normal_pdf",
    "robust_normal_params",
    "score_calls",
    "type7_quantile",
]
