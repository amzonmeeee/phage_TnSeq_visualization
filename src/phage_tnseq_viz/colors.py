"""PHROG functional-category colour scheme.

Colours follow the phold / pharokka Circos plots so that this tool's linear maps
are visually consistent with them.
Reference: https://github.com/gbouras13/phold (src/phold/plot/plot.py)

In this pre-visualization tool the category colour is used for the *edge* (border)
of each gene arrow, while the *fill* is reserved for essentiality (drawn black for
now, since no sequencing data is loaded yet).
"""

from __future__ import annotations

# Canonical PHROG category label -> hex colour.
# Keys are the lower-cased category strings as written by pharokka/phold in the
# GenBank CDS "function" qualifier.
PHROG_COLORS: dict[str, str] = {
    "unknown function": "#AAAAAA",
    "other": "#4deeea",
    "transcription regulation": "#ffe700",
    "dna, rna and nucleotide metabolism": "#f000ff",
    "lysis": "#001eff",
    "moron, auxiliary metabolic gene and host takeover": "#8900ff",
    "integration and excision": "#E0B0FF",
    "head and packaging": "#ff008d",
    "connector": "#5A5A5A",
    "tail": "#74ee15",
}

# Human-readable labels for the legend (canonical key -> pretty label).
PHROG_LABELS: dict[str, str] = {
    "unknown function": "Unknown function",
    "other": "Other function",
    "transcription regulation": "Transcription regulation",
    "dna, rna and nucleotide metabolism": "DNA/RNA & nucleotide metabolism",
    "lysis": "Lysis",
    "moron, auxiliary metabolic gene and host takeover": "Moron, auxiliary metabolic gene",
    "integration and excision": "Integration & excision",
    "head and packaging": "Head & packaging",
    "connector": "Connector",
    "tail": "Tail",
}

# Special colour for virulence factor / AMR / anti-CRISPR / defence-finder hits.
VFDB_AMR_COLOR = "#FF0000"
VFDB_AMR_LABEL = "VF / AMR / ACR / DefenseFinder"

# Edge colour / style for ORFs that were called de-novo (no annotation at all).
UNANNOTATED_EDGE = "#999999"

DEFAULT_CATEGORY = "unknown function"


def _normalize(function: str | None) -> str:
    if not function:
        return DEFAULT_CATEGORY
    key = function.strip().lower()
    if key in PHROG_COLORS:
        return key
    # A few common aliases seen across pharokka/phold versions.
    aliases = {
        "unknown": "unknown function",
        "moron, auxiliary metabolic gene and host takeover ": "moron, auxiliary metabolic gene and host takeover",
        "dna, rna and nucleotide metabolism ": "dna, rna and nucleotide metabolism",
        "other function": "other",
    }
    return aliases.get(key, DEFAULT_CATEGORY)


def category_color(function: str | None) -> str:
    """Return the border colour for a gene given its PHROG function string."""
    return PHROG_COLORS[_normalize(function)]


def category_key(function: str | None) -> str:
    """Return the canonical category key for a gene's function string."""
    return _normalize(function)


def legend_entries(categories: set[str], include_vfdb: bool = False) -> list[tuple[str, str]]:
    """Return ordered (label, colour) tuples for the categories that are present."""
    entries: list[tuple[str, str]] = []
    for key, color in PHROG_COLORS.items():
        if key in categories:
            entries.append((PHROG_LABELS[key], color))
    if include_vfdb:
        entries.append((VFDB_AMR_LABEL, VFDB_AMR_COLOR))
    return entries
