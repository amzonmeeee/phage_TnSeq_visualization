"""Canonical CSV I/O for a processed Tn-Seq final dataset.

The small :class:`FinalSite` record deliberately stays independent from the
external-processing and essentiality modules.  It is the bridge between a
user-supplied final CSV, TRANSIT TPP's parsed WIG counts, gene classification,
and the read-count plot overlay.

Canonical CSV columns are ``contig``, ``position``, and ``read_count``.  The
reader also accepts common spellings from TRANSIT/R workflows, including
``accession``/``chrom``, ``TA_site``/``insertion_site``, and ``count``/
``mean_count``.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
import math
from pathlib import Path
from typing import Iterable, Mapping, Sequence


class DatasetError(ValueError):
    """Raised when a final-dataset table is malformed."""


@dataclass(frozen=True)
class FinalSite:
    """One 1-based candidate insertion position and final read count."""

    contig: str
    position: int
    read_count: float
    raw_read_count: float | None = None
    read_count_sd: float | None = None
    n_subsamples: int | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "contig", _contig(self.contig))
        object.__setattr__(self, "position", _position(self.position))
        object.__setattr__(self, "read_count", _count(self.read_count, "read_count"))
        if self.raw_read_count is not None:
            object.__setattr__(
                self, "raw_read_count", _count(self.raw_read_count, "raw_read_count")
            )
        if self.read_count_sd is not None:
            object.__setattr__(
                self, "read_count_sd", _count(self.read_count_sd, "read_count_sd")
            )
        if self.n_subsamples is not None:
            value = _position(self.n_subsamples, "n_subsamples")
            object.__setattr__(self, "n_subsamples", value)


_CONTIG_COLUMNS = ("contig", "accession", "chrom", "chromosome")
_POSITION_COLUMNS = ("position", "insertion_site", "ta_site", "site")
_COUNT_COLUMNS = ("read_count", "count", "mean_count")
_RAW_COLUMNS = ("raw_read_count", "raw_count")
_SD_COLUMNS = ("read_count_sd", "sd_count", "count_sd")
_N_COLUMNS = ("n_subsamples", "subsamples", "n_replicates")


def load_final_dataset(
    path: str | Path,
    *,
    default_contig: str | None = None,
) -> list[FinalSite]:
    """Read, validate, and coalesce a user-supplied final-site CSV.

    A one-contig CSV may omit the contig column when ``default_contig`` is
    supplied.  Duplicate rows at one coordinate are merged by summing their
    counts, which makes ordinary count-table exports safe to import.
    """
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"final dataset not found: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise DatasetError(f"final dataset has no header: {path}")
        columns = {name.strip().lower(): name for name in reader.fieldnames if name}
        contig_column = _find_column(columns, _CONTIG_COLUMNS)
        position_column = _required_column(columns, _POSITION_COLUMNS, "position")
        count_column = _required_column(columns, _COUNT_COLUMNS, "read_count")
        raw_column = _find_column(columns, _RAW_COLUMNS)
        sd_column = _find_column(columns, _SD_COLUMNS)
        n_column = _find_column(columns, _N_COLUMNS)
        if contig_column is None and not default_contig:
            raise DatasetError(
                "final dataset needs one of contig/accession/chrom, or a default_contig"
            )

        rows: list[FinalSite] = []
        for line_number, row in enumerate(reader, start=2):
            try:
                contig = _cell(row, contig_column) if contig_column else default_contig
                rows.append(
                    FinalSite(
                        contig=contig or "",
                        position=_cell(row, position_column),
                        read_count=_cell(row, count_column),
                        raw_read_count=_optional_cell(row, raw_column),
                        read_count_sd=_optional_cell(row, sd_column),
                        n_subsamples=_optional_cell(row, n_column),
                    )
                )
            except (DatasetError, ValueError) as exc:
                raise DatasetError(f"{path}:{line_number}: {exc}") from exc
    return coalesce_final_sites(rows)


def write_final_dataset(
    path: str | Path,
    sites: Iterable[FinalSite],
    *,
    gene_assignments: Mapping[tuple[str, int], Sequence[tuple[str, str]]] | None = None,
) -> Path:
    """Write canonical final sites plus optional per-site CDS annotations.

    ``gene_assignments`` is keyed by ``(contig, position)`` and contains ordered
    ``(gene_id, strand)`` pairs.  Semicolon-delimited values preserve matching
    gene/strand order in a compact CSV; the in-memory essentiality path retains
    the fully normalised many-to-many representation.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    assignments = gene_assignments or {}
    fields = [
        "contig", "position", "raw_read_count", "read_count", "read_count_sd",
        "n_subsamples", "gene_ids", "gene_strands",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for site in coalesce_final_sites(sites):
            genes = tuple(assignments.get((site.contig, site.position), ()))
            writer.writerow(
                {
                    "contig": site.contig,
                    "position": site.position,
                    "raw_read_count": _format_optional(site.raw_read_count),
                    "read_count": _format_count(site.read_count),
                    "read_count_sd": _format_optional(site.read_count_sd),
                    "n_subsamples": "" if site.n_subsamples is None else site.n_subsamples,
                    "gene_ids": ";".join(gene_id for gene_id, _ in genes),
                    "gene_strands": ";".join(strand for _, strand in genes),
                }
            )
    return path


def final_site_from_count(count: object) -> FinalSite:
    """Convert a pipeline ``InsertionCount`` or averaged-count-like object."""
    contig = getattr(count, "contig")
    position = getattr(count, "position")
    if hasattr(count, "mean_count"):
        read_count = getattr(count, "mean_count")
        raw_read_count = None
        read_count_sd = getattr(count, "sd_count", None)
        n_subsamples = getattr(count, "n_subsamples", None)
    else:
        read_count = getattr(count, "read_count", None)
        if read_count is None:
            read_count = getattr(count, "count")
        raw_read_count = read_count
        read_count_sd = None
        n_subsamples = None
    return FinalSite(
        contig=contig,
        position=position,
        read_count=read_count,
        raw_read_count=raw_read_count,
        read_count_sd=read_count_sd,
        n_subsamples=n_subsamples,
    )


def final_sites_from_counts(counts: Iterable[object]) -> list[FinalSite]:
    """Convert and coalesce parsed TPP counts or averaged subsample counts."""
    return coalesce_final_sites(final_site_from_count(count) for count in counts)


def coalesce_final_sites(sites: Iterable[FinalSite]) -> list[FinalSite]:
    """Merge duplicate coordinate rows by summing counts deterministically."""
    grouped: dict[tuple[str, int], list[FinalSite]] = {}
    for site in sites:
        if not isinstance(site, FinalSite):
            site = FinalSite(**site)  # type: ignore[arg-type]
        grouped.setdefault((site.contig, site.position), []).append(site)

    merged: list[FinalSite] = []
    for (contig, position), rows in sorted(grouped.items()):
        raw_values = [row.raw_read_count for row in rows if row.raw_read_count is not None]
        sd_values = [row.read_count_sd for row in rows if row.read_count_sd is not None]
        n_values = {row.n_subsamples for row in rows if row.n_subsamples is not None}
        merged.append(
            FinalSite(
                contig=contig,
                position=position,
                read_count=sum(row.read_count for row in rows),
                raw_read_count=sum(raw_values) if raw_values else None,
                # Independent duplicate records are summed, so their standard
                # deviations combine in quadrature.  This is only metadata; the
                # final count remains the source of truth for classification.
                read_count_sd=math.sqrt(sum(value * value for value in sd_values)) if sd_values else None,
                n_subsamples=n_values.pop() if len(n_values) == 1 else None,
            )
        )
    return merged


def fill_missing_final_sites(
    sites: Iterable[FinalSite],
    candidate_sites: Mapping[str, Iterable[int]],
) -> list[FinalSite]:
    """Add explicit zero-count records for every omitted candidate coordinate."""
    existing = {(site.contig, site.position): site for site in coalesce_final_sites(sites)}
    for contig, positions in candidate_sites.items():
        for position in positions:
            key = (_contig(contig), _position(position))
            existing.setdefault(key, FinalSite(key[0], key[1], 0.0, raw_read_count=0.0))
    return [existing[key] for key in sorted(existing)]


def group_counts_for_plotting(sites: Iterable[FinalSite]) -> dict[str, dict[int, float]]:
    """Return the ``accession -> position -> final count`` shape used by render."""
    return {
        contig: {position: site.read_count for (site_contig, position), site in grouped.items()
                 if site_contig == contig}
        for contig, grouped in _group_by_contig(coalesce_final_sites(sites)).items()
    }


def _group_by_contig(sites: Iterable[FinalSite]) -> dict[str, dict[tuple[str, int], FinalSite]]:
    grouped: dict[str, dict[tuple[str, int], FinalSite]] = {}
    for site in sites:
        grouped.setdefault(site.contig, {})[(site.contig, site.position)] = site
    return grouped


def _find_column(columns: Mapping[str, str], aliases: Sequence[str]) -> str | None:
    return next((columns[alias] for alias in aliases if alias in columns), None)


def _required_column(columns: Mapping[str, str], aliases: Sequence[str], label: str) -> str:
    value = _find_column(columns, aliases)
    if value is None:
        raise DatasetError(f"final dataset needs a {label} column (accepted: {', '.join(aliases)})")
    return value


def _cell(row: Mapping[str, str | None], column: str | None) -> str:
    if column is None:
        return ""
    value = row.get(column)
    if value is None or not str(value).strip():
        raise DatasetError(f"missing value for {column!r}")
    return str(value).strip()


def _optional_cell(row: Mapping[str, str | None], column: str | None) -> str | None:
    if column is None:
        return None
    value = row.get(column)
    return str(value).strip() if value is not None and str(value).strip() else None


def _contig(value: object) -> str:
    result = str(value).strip() if value is not None else ""
    if not result:
        raise DatasetError("contig must not be blank")
    return result


def _position(value: object, label: str = "position") -> int:
    if isinstance(value, bool):
        raise DatasetError(f"{label} must be a positive integer")
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise DatasetError(f"{label} must be a positive integer") from exc
    if not math.isfinite(number) or not number.is_integer() or number < 1:
        raise DatasetError(f"{label} must be a positive integer")
    return int(number)


def _count(value: object, label: str) -> float:
    if isinstance(value, bool):
        raise DatasetError(f"{label} must be a non-negative number")
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise DatasetError(f"{label} must be a non-negative number") from exc
    if not math.isfinite(number) or number < 0:
        raise DatasetError(f"{label} must be a non-negative number")
    return number


def _format_count(value: float) -> str:
    return f"{value:g}"


def _format_optional(value: float | None) -> str:
    return "" if value is None else _format_count(value)


__all__ = [
    "DatasetError",
    "FinalSite",
    "coalesce_final_sites",
    "fill_missing_final_sites",
    "final_site_from_count",
    "final_sites_from_counts",
    "group_counts_for_plotting",
    "load_final_dataset",
    "write_final_dataset",
]
