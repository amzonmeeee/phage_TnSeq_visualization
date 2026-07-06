"""Optional genomic feature tracks: GC content, GC skew, and non-CDS features.

GC content / skew are computed in sliding windows following the same conventions
used by pharokka/phold plots (deviation from the genome-wide GC fraction).
tRNA / tmRNA / CRISPR features are read straight from the GenBank record.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from Bio.SeqRecord import SeqRecord


@dataclass
class WindowTrack:
    """A per-window numeric signal (GC content or GC skew)."""

    positions: list[float] = field(default_factory=list)  # window mid-points (1-based)
    values: list[float] = field(default_factory=list)     # deviation from genome mean


@dataclass
class NonCdsFeature:
    start: int
    end: int
    strand: int
    kind: str          # "tRNA", "tmRNA", "CRISPR", ...
    label: str | None = None


def _gc_fraction(seq: str) -> float:
    if not seq:
        return 0.0
    g = seq.count("G") + seq.count("g")
    c = seq.count("C") + seq.count("c")
    return (g + c) / len(seq)


def compute_gc_content(sequence: str, window: int, step: int) -> WindowTrack:
    """GC content per window, expressed as deviation from the genome-wide GC."""
    seq = sequence.upper()
    n = len(seq)
    if n == 0:
        return WindowTrack()
    genome_gc = _gc_fraction(seq)
    track = WindowTrack()
    half = window // 2
    for start in range(0, n, step):
        sub = seq[start : start + window]
        if not sub:
            continue
        track.positions.append(start + half + 1)
        track.values.append(_gc_fraction(sub) - genome_gc)
    return track


def compute_gc_skew(sequence: str, window: int, step: int) -> WindowTrack:
    """GC skew (G-C)/(G+C) per window."""
    seq = sequence.upper()
    n = len(seq)
    if n == 0:
        return WindowTrack()
    track = WindowTrack()
    half = window // 2
    for start in range(0, n, step):
        sub = seq[start : start + window]
        if not sub:
            continue
        g = sub.count("G")
        c = sub.count("C")
        skew = (g - c) / (g + c) if (g + c) else 0.0
        track.positions.append(start + half + 1)
        track.values.append(skew)
    return track


def compute_insertion_density(
    sites: list[int], length: int, window: int, step: int | None = None
) -> WindowTrack:
    """Insertion-site density in insertions-per-kb over sliding windows.

    ``positions`` are window mid-points (1-based); ``values`` are sites/kb.
    """
    if length <= 0:
        return WindowTrack()
    from bisect import bisect_right

    step = step or max(1, window // 2)
    sites = sorted(sites)
    track = WindowTrack()
    half = window // 2
    for start in range(0, length, step):
        end = min(start + window, length)
        w = end - start
        if w <= 0:
            continue
        # count sites in (start, end]  (1-based sites)
        count = bisect_right(sites, end) - bisect_right(sites, start)
        track.positions.append(start + half + 1)
        track.values.append(count / w * 1000.0)  # per kb
    return track


_NONCDS_TYPES = {"tRNA": "tRNA", "tmRNA": "tmRNA", "ncRNA": "ncRNA",
                 "CRISPR": "CRISPR", "repeat_region": "CRISPR"}


def extract_noncds_features(record: SeqRecord) -> list[NonCdsFeature]:
    """Pull tRNA / tmRNA / CRISPR features from a GenBank record."""
    feats: list[NonCdsFeature] = []
    for feat in record.features:
        kind = _NONCDS_TYPES.get(feat.type)
        if kind is None:
            continue
        # Only treat repeat_region as CRISPR when it looks like one.
        if feat.type == "repeat_region":
            note = " ".join(feat.qualifiers.get("note", []) + feat.qualifiers.get("rpt_family", []))
            if "crispr" not in note.lower():
                continue
        start = int(feat.location.start) + 1
        end = int(feat.location.end)
        strand = 1 if (feat.location.strand or 1) >= 0 else -1
        label = None
        for k in ("product", "note", "gene"):
            if feat.qualifiers.get(k):
                label = str(feat.qualifiers[k][0])
                break
        feats.append(NonCdsFeature(start, end, strand, kind, label))
    return feats
