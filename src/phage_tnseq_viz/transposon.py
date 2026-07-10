"""Transposon insertion-site preferences and motif scanning.

A transposon's insertion preference is expressed as a short DNA motif using IUPAC
degenerate bases (e.g. mariner/Himar1 inserts at ``TA``). This module turns such a
motif into genomic positions where the transposon can insert.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

# Built-in presets: name -> IUPAC insertion motif.
# ``None`` means the transposon has little/no sequence preference (≈ random), so we
# do not draw insertion ticks for it (that would just fill the whole genome).
PRESETS: dict[str, str | None] = {
    "mariner": "TA",   # Himar1 / mariner family strictly insert at TA dinucleotides
    "himar1": "TA",
    "tn5": None,        # very weak preference, treated as random here
    "tn10": None,       # weak "NGCTNAGCN"-like consensus; treat as random by default
    "tn7": None,        # site-specific (attTn7), not motif-based
    "mu": None,         # near-random
}

# IUPAC degenerate base -> regex character class.
_IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]",
    "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]",
    "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]",
}

_COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")


@dataclass
class Transposon:
    name: str
    motif: str | None  # IUPAC motif, or None for random / site-specific

    @property
    def has_preference(self) -> bool:
        return bool(self.motif)


def resolve_transposon(spec: str) -> Transposon:
    """Resolve a --transposon argument.

    ``spec`` may be a preset name (e.g. ``mariner``) or a literal IUPAC motif
    (e.g. ``TA`` or ``NTAN``). Presets take precedence.
    """
    key = spec.strip().lower()
    if key in PRESETS:
        return Transposon(name=key, motif=PRESETS[key])
    motif = spec.strip().upper()
    _validate_motif(motif)
    return Transposon(name=f"custom:{motif}", motif=motif)


def _validate_motif(motif: str) -> None:
    if not motif:
        raise ValueError("Empty transposon motif.")
    bad = sorted(set(motif) - set(_IUPAC))
    if bad:
        raise ValueError(
            f"Invalid IUPAC base(s) in motif {motif!r}: {''.join(bad)}. "
            f"Allowed: {''.join(_IUPAC)}"
        )


def _motif_to_regex(motif: str) -> str:
    return "".join(_IUPAC[b] for b in motif)


def _reverse_complement(motif: str) -> str:
    return motif.translate(_COMPLEMENT)[::-1]


def find_insertion_sites(
    sequence: str,
    motif: str,
    *,
    both_strands: bool = True,
) -> list[int]:
    """Return sorted 1-based genomic positions of the motif's insertion point.

    The insertion point is the motif's 5' (first) base *on the strand it occurs on*.
    Overlapping matches are found. If ``both_strands`` is set the reverse-complement
    motif is also scanned on the minus strand, and its 5' base maps to the *right*
    end of the match in forward coordinates (so the point stays on the same base of
    the motif regardless of strand). Hits are deduplicated across strands; for a
    palindromic motif such as ``TA`` the minus-strand scan is skipped entirely
    because it would only re-find the same sites.
    """
    seq = sequence.upper()
    n = len(motif)
    positions: set[int] = set()

    fwd = re.compile(f"(?=({_motif_to_regex(motif)}))")
    for m in fwd.finditer(seq):
        positions.add(m.start() + 1)  # 1-based; 5' base of the motif on the + strand

    if both_strands:
        rc = _reverse_complement(motif)
        if rc != motif:
            rev = re.compile(f"(?=({_motif_to_regex(rc)}))")
            for m in rev.finditer(seq):
                # The rev-comp match spans forward 1-based positions
                # [m.start()+1 .. m.start()+n]. The motif's 5' base sits on the minus
                # strand, i.e. at the right end in forward coordinates.
                positions.add(m.start() + n)

    return sorted(positions)
