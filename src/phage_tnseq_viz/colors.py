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

# Gene-arrow fill colours for the built-in Harms-lab saturation classifier.
#
# PHROG stays the *border* colour; the arrow *fill* communicates the Tn-Seq call.
# The biological calls run as a single warm sequential ramp from deep red
# (essential) through orange and yellow to white (non-essential) — depletion reads
# as "hotter", tolerance fades to blank — so no out-of-scale hue (the old green /
# blue) interrupts the ramp.  The non-call states stay a neutral grey, kept off the
# ramp so "no data" is never mistaken for a fitness level.
ESSENTIALITY_COLORS: dict[str, str] = {
    "Essential": "#d73027",
    "Strong fitness defect": "#fc8d59",
    "Intermediate": "#fee08b",
    "Reduced fitness": "#ffffbf",
    "Non-essential": "#fffde0",
    "Ambiguous": "#9e9e9e",
    "Insufficient sites": "#d9d9d9",
    "Unclassified": "#bdbdbd",
}

ESSENTIALITY_LABELS: dict[str, str] = {
    "Essential": "Essential",
    "Strong fitness defect": "Strong fitness defect",
    "Intermediate": "Intermediate",
    "Reduced fitness": "Reduced fitness",
    "Non-essential": "Non-essential",
    "Ambiguous": "Ambiguous",
    "Insufficient sites": "Insufficient candidate sites (≤5)",
    "Unclassified": "Unclassified / custom call",
}


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


# Ordered keyword rules to recover a PHROG category from a plain product string.
# Needed for GenBank files that carry only `/product` (e.g. NCBI downloads, or a
# phold/phynteny file round-tripped through Benchling, which drops `/function`).
# Order matters: the first category whose keyword is found in the product wins,
# so more specific structural terms are listed before generic metabolism ones.
_PRODUCT_RULES: list[tuple[str, tuple[str, ...]]] = [
    ("integration and excision",
     ("integrase", "recombinase", "excisionase", "transposase", "resolvase",
      "site-specific", "invertase")),
    ("lysis",
     ("holin", "endolysin", "spanin", "lysin", "lysozyme", "amidase",
      "peptidoglycan", "murein", "cell wall", "lysis")),
    # Connector is checked before tail so "head-tail adaptor" is not caught by
    # the "tail" keyword.
    ("connector",
     ("connector", "adaptor", "adapter", "head-tail", "head to tail", "neck",
      "collar", "stopper")),
    ("tail",
     ("tail", "fiber", "fibre", "baseplate", "base plate", "tape measure",
      "tape-measure", "tapemeasure", "sheath", "spike", "receptor binding",
      "receptor-binding")),
    ("head and packaging",
     ("terminase", "capsid", "portal", "prohead", "scaffold", "head", "packaging",
      "coat protein", "maturation protease")),
    ("transcription regulation",
     ("transcriptional regulator", "transcription regulator", "repressor",
      "anti-repressor", "antirepressor", "sigma factor", "rna polymerase",
      "transcription")),
    ("dna, rna and nucleotide metabolism",
     ("polymerase", "primase", "helicase", "exonuclease", "endonuclease",
      "nuclease", "ligase", "kinase", "methylase", "methyltransferase",
      "phosphatase", "reductase", "ribonucleotide", "thymidylate", "dutpase",
      "topoisomerase", "single-strand", "single strand", "recombination",
      "replication", "nucleotide", "rnr", "dna ", "rna ", "deaminase",
      "glycosyltransferase")),
    ("moron, auxiliary metabolic gene and host takeover",
     ("antitoxin", "toxin", "moron", "host takeover", "superinfection",
      "immunity", "membrane protein")),
]

_UNKNOWN_KEYS = ("hypothetical", "unknown", "uncharacteri", "putative protein",
                 "duf", "domain of unknown")


def infer_category_from_product(product: str | None) -> str | None:
    """Best-effort PHROG category from a free-text product name.

    Returns a canonical category key (a key of ``PHROG_COLORS``) or ``None`` when
    nothing matches confidently (caller then falls back to "unknown function").
    """
    if not product:
        return None
    text = product.strip().lower()
    if not text:
        return None
    for key in _UNKNOWN_KEYS:
        if key in text:
            return "unknown function"
    for category, keywords in _PRODUCT_RULES:
        for kw in keywords:
            if kw in text:
                return category
    return None


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


def essentiality_key(call: str | None) -> str:
    """Normalise a built-in essentiality call for plotting.

    The original R script writes ``NA`` for genes with five or fewer candidate
    sites.  Python represents that as ``None`` and the plot labels it explicitly
    as ``"Insufficient sites"``.  A caller-supplied classifier may return any
    other label; those are displayed as ``Unclassified`` unless it supplies its
    own palette in a future extension.
    """
    if call is None or not str(call).strip():
        return "Insufficient sites"
    value = str(call).strip()
    return value if value in ESSENTIALITY_COLORS else "Unclassified"


def essentiality_color(call: str | None) -> str:
    """Return the arrow-fill colour for an essentiality call."""
    return ESSENTIALITY_COLORS[essentiality_key(call)]


def essentiality_label(call: str | None) -> str:
    """Return the legend label for an essentiality call."""
    return ESSENTIALITY_LABELS[essentiality_key(call)]


def essentiality_legend_entries(calls: set[str | None]) -> list[tuple[str, str]]:
    """Return present essentiality calls in a stable biological order."""
    present = {essentiality_key(call) for call in calls}
    return [
        (ESSENTIALITY_LABELS[key], colour)
        for key, colour in ESSENTIALITY_COLORS.items()
        if key in present
    ]
