"""Library quality metrics for a per-site count table.

Every essentiality call downstream rests on assumptions about the library:
that enough candidate sites were hit to distinguish an empty gene from an
unlucky one, and that the counts are not dominated by a handful of runaway
sites.  These metrics report on those assumptions directly, so a result can be
judged before it is believed.

The metric set follows TRANSIT's ``tnseq_stats``, which is the established
vocabulary for Tn-Seq library QC, and the interpretation thresholds in
:data:`WARNING_RULES` are its published guidance.  Two things differ, both
because this tool targets phage rather than bacterial genomes:

* Metrics are reported per contig, matching how saturation and read thresholds
  are already computed elsewhere in this package.
* The Pickands tail index is adapted to small genomes; see
  :func:`pickands_tail_index`.

Everything here is pure Python, like the rest of the analysis modules, so the
QC path adds no dependency.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass
import math


# TRANSIT's published guidance for spotting a problem dataset.  These were
# derived from bacterial libraries; on a phage genome, with far fewer candidate
# sites, treat them as prompts to look rather than as pass/fail limits.
WARNING_RULES: tuple[tuple[str, str], ...] = (
    ("density", "saturation below 0.30: many genes will be uncallable"),
    ("nz_mean", "mean count at hit sites below 10: library may be under-sequenced"),
    ("max_count", "one site exceeds 1e6 reads: likely an outlier distorting the scale"),
    ("skewness", "skewness above 50: count distribution is very noisy"),
    ("pickands_tail_index", "Pickands tail index above 0.5: heavy-tailed counts"),
)

_DENSITY_FLOOR = 0.30
_NZ_MEAN_FLOOR = 10.0
_MAX_COUNT_CEILING = 1e6
_SKEWNESS_CEILING = 50.0
_PTI_CEILING = 0.5

# Pickands compares order statistics at ranks M, 2M and 4M.  TRANSIT scans
# M = 10..99, which silently requires at least 397 sites -- more than a small
# phage genome has.  The scan is therefore capped by the data instead, which
# reproduces TRANSIT exactly whenever the data is large enough for its range.
_PTI_MIN_M = 10
_PTI_MAX_M = 100


@dataclass(frozen=True)
class DatasetStats:
    """Quality metrics for one contig's candidate-site counts.

    ``skewness``, ``kurtosis`` and ``pickands_tail_index`` are ``None`` when the
    data cannot support them: the moments need at least two hit sites with some
    spread, and the tail index needs enough sites to compare order statistics at
    ranks M, 2M and 4M.
    """

    dataset: str
    sites: int
    hit_sites: int
    density: float
    mean_count: float
    nz_mean: float
    nz_median: float
    max_count: float
    total_counts: float
    skewness: float | None
    kurtosis: float | None
    pickands_tail_index: float | None

    @property
    def warnings(self) -> tuple[str, ...]:
        """Published thresholds this contig trips, in reporting order."""

        found: list[str] = []
        if self.density < _DENSITY_FLOOR:
            found.append(f"saturation {self.density:.2f} is below {_DENSITY_FLOOR:.2f}")
        if self.hit_sites and self.nz_mean < _NZ_MEAN_FLOOR:
            found.append(f"mean count at hit sites {self.nz_mean:.1f} is below {_NZ_MEAN_FLOOR:.0f}")
        if self.max_count > _MAX_COUNT_CEILING:
            found.append(f"top site has {self.max_count:.0f} reads, a likely outlier")
        if self.skewness is not None and self.skewness > _SKEWNESS_CEILING:
            found.append(f"skewness {self.skewness:.1f} is above {_SKEWNESS_CEILING:.0f}")
        if self.pickands_tail_index is not None and self.pickands_tail_index > _PTI_CEILING:
            found.append(
                f"Pickands tail index {self.pickands_tail_index:.2f} is above {_PTI_CEILING:.1f}"
            )
        return tuple(found)


def central_moment_skewness(values: Sequence[float]) -> float | None:
    """Sample skewness ``m3 / m2**1.5`` (the biased ``g1`` estimator).

    Returns ``None`` for fewer than two values or for a constant sequence, where
    skewness is undefined rather than zero.
    """

    n = len(values)
    if n < 2:
        return None
    mean = sum(values) / n
    m2 = sum((value - mean) ** 2 for value in values) / n
    if m2 <= 0.0:
        return None
    m3 = sum((value - mean) ** 3 for value in values) / n
    return m3 / m2**1.5


def central_moment_kurtosis(values: Sequence[float]) -> float | None:
    """Sample *excess* kurtosis ``m4 / m2**2 - 3``, so a normal sample gives 0."""

    n = len(values)
    if n < 2:
        return None
    mean = sum(values) / n
    m2 = sum((value - mean) ** 2 for value in values) / n
    if m2 <= 0.0:
        return None
    m4 = sum((value - mean) ** 4 for value in values) / n
    return m4 / m2**2 - 3.0


def pickands_tail_index(counts: Sequence[float]) -> float | None:
    """Median Pickands estimate of the tail index; higher means heavier tails.

    For counts sorted descending, each estimate is
    ``log((X[M] - X[2M]) / (X[2M] - X[4M])) / log(2)``, and the median over a
    range of ``M`` is reported.  A heavier tail means a few sites carry a
    disproportionate share of the reads, which is what makes a read-count scale
    unreadable and what ``--read-histogram-cap`` exists to mitigate.

    TRANSIT scans ``M = 10..99``, which needs at least 397 sites.  Phage genomes
    routinely have fewer, so the scan is capped at ``(n - 1) // 4`` here.  When
    the data is large enough for TRANSIT's full range the two agree exactly;
    when fewer than 41 sites are available no estimate is possible and ``None``
    is returned.

    Estimates where the order statistics tie -- common once the tail runs into
    the zero counts -- are undefined and are skipped rather than poisoning the
    median.
    """

    n = len(counts)
    highest_m = min(_PTI_MAX_M - 1, (n - 1) // 4)
    if highest_m < _PTI_MIN_M:
        return None

    ordered = sorted(counts, reverse=True)
    estimates: list[float] = []
    for m in range(_PTI_MIN_M, highest_m + 1):
        upper = ordered[m] - ordered[2 * m]
        lower = ordered[2 * m] - ordered[4 * m]
        if upper <= 0.0 or lower <= 0.0:
            continue
        estimates.append(math.log(upper / lower) / math.log(2.0))
    if not estimates:
        return None
    return _median(estimates)


def dataset_stats(counts: Iterable[float], *, dataset: str = "") -> DatasetStats:
    """Compute every metric for one contig's complete candidate-site counts.

    ``counts`` must include the zero-count sites; without them saturation and
    the mean over all sites are both meaningless.
    """

    values = [float(count) for count in counts]
    n = len(values)
    if n == 0:
        return DatasetStats(
            dataset=dataset, sites=0, hit_sites=0, density=0.0, mean_count=0.0,
            nz_mean=0.0, nz_median=0.0, max_count=0.0, total_counts=0.0,
            skewness=None, kurtosis=None, pickands_tail_index=None,
        )

    positive = [value for value in values if value > 0.0]
    total = sum(values)
    return DatasetStats(
        dataset=dataset,
        sites=n,
        hit_sites=len(positive),
        density=len(positive) / n,
        mean_count=total / n,
        nz_mean=(sum(positive) / len(positive)) if positive else 0.0,
        nz_median=_median(positive) if positive else 0.0,
        max_count=max(values),
        total_counts=total,
        # Moments are taken over hit sites only: the zeros are already reported
        # by density, and including them would make every sparse library look
        # skewed for a reason that has nothing to do with count noise.
        skewness=central_moment_skewness(positive),
        kurtosis=central_moment_kurtosis(positive),
        pickands_tail_index=pickands_tail_index(values),
    )


def stats_by_contig(
    rows: Iterable[tuple[str, float]],
) -> list[DatasetStats]:
    """Compute per-contig metrics from ``(contig, read_count)`` pairs."""

    grouped: dict[str, list[float]] = {}
    for contig, count in rows:
        grouped.setdefault(contig, []).append(float(count))
    return [dataset_stats(counts, dataset=contig) for contig, counts in sorted(grouped.items())]


TABLE_COLUMNS: tuple[str, ...] = (
    "dataset", "sites", "hit_sites", "density", "mean_count", "nz_mean",
    "nz_median", "max_count", "total_counts", "skewness", "kurtosis",
    "pickands_tail_index",
)


def as_row(stats: DatasetStats) -> dict[str, str]:
    """Render one result as CSV-ready strings, blank where a metric is undefined."""

    def number(value: float | None, spec: str) -> str:
        return "" if value is None else format(value, spec)

    return {
        "dataset": stats.dataset,
        "sites": str(stats.sites),
        "hit_sites": str(stats.hit_sites),
        "density": f"{stats.density:.3f}",
        "mean_count": f"{stats.mean_count:.1f}",
        "nz_mean": f"{stats.nz_mean:.1f}",
        "nz_median": f"{stats.nz_median:.1f}",
        "max_count": f"{stats.max_count:.0f}",
        "total_counts": f"{stats.total_counts:.0f}",
        "skewness": number(stats.skewness, ".2f"),
        "kurtosis": number(stats.kurtosis, ".2f"),
        "pickands_tail_index": number(stats.pickands_tail_index, ".3f"),
    }


def format_table(results: Sequence[DatasetStats]) -> list[str]:
    """Render an aligned plain-text table, one line per element."""

    if not results:
        return []
    headers = ("contig", "sites", "hit", "density", "NZmean", "NZmedian", "max", "skew", "PTI")
    rows = [
        (
            stats.dataset or "-",
            str(stats.sites),
            str(stats.hit_sites),
            f"{stats.density:.3f}",
            f"{stats.nz_mean:.1f}",
            f"{stats.nz_median:.1f}",
            f"{stats.max_count:.0f}",
            "-" if stats.skewness is None else f"{stats.skewness:.1f}",
            "-" if stats.pickands_tail_index is None else f"{stats.pickands_tail_index:.2f}",
        )
        for stats in results
    ]
    widths = [max(len(header), *(len(row[i]) for row in rows)) for i, header in enumerate(headers)]
    lines = ["  ".join(header.rjust(widths[i]) for i, header in enumerate(headers))]
    lines.extend("  ".join(cell.rjust(widths[i]) for i, cell in enumerate(row)) for row in rows)
    return lines


def _median(values: Sequence[float]) -> float:
    ordered = sorted(values)
    n = len(ordered)
    middle = n // 2
    if n % 2:
        return ordered[middle]
    return (ordered[middle - 1] + ordered[middle]) / 2.0


__all__ = [
    "TABLE_COLUMNS",
    "WARNING_RULES",
    "DatasetStats",
    "as_row",
    "central_moment_kurtosis",
    "central_moment_skewness",
    "dataset_stats",
    "format_table",
    "pickands_tail_index",
    "stats_by_contig",
]
