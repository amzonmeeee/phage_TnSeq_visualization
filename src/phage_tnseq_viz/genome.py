"""Load a phage genome from GenBank and extract genes.

Two input flavours are supported:

1. Annotated GenBank produced by pharokka / phold / phynteny(_transformer): CDS
   features carry a ``function`` qualifier (PHROG category) and ``product``.
2. A plain GenBank download with no CDS annotation: ORFs are called de-novo with
   pyrodigal-gv (the same viral gene caller pharokka uses internally).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from .colors import infer_category_from_product
from .tracks import NonCdsFeature, extract_noncds_features


@dataclass
class Gene:
    start: int          # 1-based inclusive
    end: int            # 1-based inclusive
    strand: int         # +1 / -1
    annotated: bool     # True if from GenBank CDS, False if called de-novo
    function: str | None = None   # PHROG category string (annotated only)
    product: str | None = None
    locus: str | None = None
    is_vfdb_amr: bool = False     # virulence/AMR/anti-CRISPR/defence hit


@dataclass
class GenomeRecord:
    name: str
    accession: str
    length: int
    sequence: str
    genes: list[Gene] = field(default_factory=list)
    noncds: list[NonCdsFeature] = field(default_factory=list)
    annotation_source: str = "genbank"  # "genbank" or "pyrodigal-gv"

    @property
    def n_genes(self) -> int:
        return len(self.genes)


def gene_identifier(gene: Gene, *, contig: str | None = None) -> str:
    """Return the stable identifier used to join genes to Tn-Seq results.

    GenBank ``locus_tag`` (or the ``ID`` / ``label`` already stored in
    :attr:`Gene.locus`) is preferred.  Some phage annotations do not supply one;
    a coordinate-and-strand fallback is deterministic within a contig and keeps
    those CDS features usable in a final-dataset CSV.  Pass ``contig`` when an
    identifier might be stored outside a per-contig mapping, as the essentiality
    module does.
    """
    if gene.locus:
        return gene.locus
    prefix = f"{contig}:" if contig else ""
    return f"{prefix}{gene.start}-{gene.end}:{gene.strand}"


# Qualifier keys that, when present, flag a CDS as a VF/AMR/ACR/DefenseFinder hit.
_VFDB_AMR_KEYS = ("vfdb", "card", "amr", "acr", "defensefinder", "acrdb")


def _get_qual(feature, key: str) -> str | None:
    vals = feature.qualifiers.get(key)
    if vals:
        return str(vals[0])
    return None


def _extract_annotated_genes(record: SeqRecord) -> list[Gene]:
    genes: list[Gene] = []
    for feat in record.features:
        if feat.type != "CDS":
            continue
        start = int(feat.location.start) + 1  # 0-based -> 1-based
        end = int(feat.location.end)
        strand = 1 if (feat.location.strand or 1) >= 0 else -1
        function = _get_qual(feat, "function")
        product = _get_qual(feat, "product")
        # Files that carry only /product (NCBI downloads, or a phold/phynteny
        # file round-tripped through Benchling, which strips /function): recover
        # the PHROG category from the product text so arrows are still coloured.
        if function is None:
            function = infer_category_from_product(product)
        locus = (_get_qual(feat, "locus_tag") or _get_qual(feat, "ID")
                 or _get_qual(feat, "label"))
        is_hit = any(
            k in feat.qualifiers or _get_qual(feat, k) for k in _VFDB_AMR_KEYS
        )
        genes.append(
            Gene(
                start=start,
                end=end,
                strand=strand,
                annotated=True,
                function=function,
                product=product,
                locus=locus,
                is_vfdb_amr=is_hit,
            )
        )
    return genes


def _call_orfs(sequence: str) -> list[Gene]:
    """Call ORFs de-novo with pyrodigal-gv."""
    try:
        import pyrodigal_gv
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError(
            "This genome has no CDS annotation and pyrodigal-gv is not installed. "
            "Install it with `pip install pyrodigal-gv`, or supply an annotated GenBank."
        ) from exc

    finder = pyrodigal_gv.ViralGeneFinder(meta=True)
    genes: list[Gene] = []
    for pred in finder.find_genes(sequence.encode()):
        genes.append(
            Gene(
                start=pred.begin,      # pyrodigal is already 1-based inclusive
                end=pred.end,
                strand=1 if pred.strand >= 0 else -1,
                annotated=False,
                function=None,
                product="hypothetical protein (de-novo ORF)",
            )
        )
    return genes


def load_genome(
    path: str | Path,
    *,
    name_override: str | None = None,
    call_orfs_if_missing: bool = True,
) -> list[GenomeRecord]:
    """Parse a GenBank file into one or more :class:`GenomeRecord`.

    One :class:`GenomeRecord` is produced per GenBank record (contig).
    """
    path = Path(path)
    records: list[GenomeRecord] = []
    for rec in SeqIO.parse(str(path), "genbank"):
        seq = str(rec.seq)
        genes = _extract_annotated_genes(rec)
        source = "genbank"
        if not genes and call_orfs_if_missing:
            genes = _call_orfs(seq)
            source = "pyrodigal-gv"

        accession = rec.id or rec.name or "unknown"
        display_name = name_override or (rec.description or rec.name or accession)
        records.append(
            GenomeRecord(
                name=display_name,
                accession=accession,
                length=len(seq),
                sequence=seq,
                genes=genes,
                noncds=extract_noncds_features(rec),
                annotation_source=source,
            )
        )

    if not records:
        raise ValueError(f"No GenBank records found in {path}")
    return records
