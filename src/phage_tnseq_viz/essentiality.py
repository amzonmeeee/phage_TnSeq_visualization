"""Essentiality calls for per-site phage Tn-Seq counts.

This module is a dependency-free port of the supplied
``essentiality_classification.R`` workflow.  It deliberately works with a
*complete* candidate-site table: an unobserved candidate insertion must be
present with ``read_count=0`` or saturation cannot be calculated correctly.

The reference R script calls those candidates ``TA_site``.  The neutral name
``position`` is used here so a caller can explicitly choose another candidate
site model (for example, every base for a near-random transposon).  That is an
adaptation of the TA-site method, not a claim that its biological thresholds
have been validated for every transposon.

Custom classifiers
------------------
Users may supply a local Python file containing exactly this callable::

    def classify_gene(gene_id, site_rows):
        # site_rows is a tuple of GeneSiteRow objects, including zero-count sites
        return "My category"  # or None for no call

Load it with :func:`load_classifier` and apply it with
:func:`apply_classifier`.  Loading a classifier executes user-provided Python,
so only use a file from a trusted source.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass, field, replace
import hashlib
import importlib.util
import math
from pathlib import Path
from typing import Any, Optional

from . import gaps


# Labels emitted by the supplied R implementation.
ESSENTIAL = "Essential"
NON_ESSENTIAL = "Non-essential"
INTERMEDIATE = "Intermediate"
STRONG_FITNESS_DEFECT = "Strong fitness defect"
REDUCED_FITNESS = "Reduced fitness"
AMBIGUOUS = "Ambiguous"
FURTHER_ANALYSIS = "Further analysis"

# Supplement to the R rules, not part of them: see binomial_min_sites().
ESSENTIAL_BINOMIAL = "Essential (binomial)"

# Significance level for the binomial short-gene supplement.
BINOMIAL_ALPHA = 0.05

# Significance level for the gap analysis, applied to BH-adjusted p-values.
GAP_ALPHA = 0.05

# Calls that already say "this gene cannot tolerate insertion".
_ESSENTIAL_CALLS = frozenset({ESSENTIAL, ESSENTIAL_BINOMIAL})


def binomial_min_sites(saturation: float, *, alpha: float = BINOMIAL_ALPHA) -> float | None:
    """Fewest candidate sites at which "no insertions at all" is significant.

    The R rules give no call to a gene with five or fewer candidate sites, which
    on a phage genome silences a large share of the annotation: genes are short,
    so TA sites per gene are few.  Yet a gene with *zero* insertions is the
    clearest possible essentiality signal, and whether that is surprising depends
    on how saturated the library is, not on an arbitrary site count.

    Under a Bernoulli model with a library-wide non-insertion probability
    ``phi = 1 - saturation``, a gene of ``n`` candidate sites is empty by chance
    with probability ``phi**n``.  Requiring ``phi**n < alpha`` gives the returned
    threshold ``log(alpha) / log(phi)``; a gene with more sites than that, and no
    insertions, is essential at that significance level.

    Returns ``None`` when no threshold can exist: a fully saturated library
    (``phi <= 0``) leaves no empty genes to call, and a library with no
    insertions at all (``phi >= 1``) makes every gene look empty.

    This follows the supplement Choudhery et al. (2021) added to TRANSIT's Gumbel
    method, where such genes are reported as ``EB`` alongside Gumbel's ``E``.
    """

    if not 0.0 <= saturation <= 1.0:
        raise ValueError(f"saturation must be between 0 and 1: {saturation!r}")
    if not 0.0 < alpha < 1.0:
        raise ValueError(f"alpha must be between 0 and 1 (exclusive): {alpha!r}")
    phi = 1.0 - saturation
    if phi <= 0.0 or phi >= 1.0:
        return None
    return math.log(alpha) / math.log(phi)


@dataclass(frozen=True)
class InsertionSite:
    """A candidate genomic insertion site and its post-processing read count.

    ``position`` is 1-based and inclusive, matching :class:`~.genome.Gene`.
    ``contig`` should normally be the GenBank record accession.
    """

    position: int
    read_count: float
    contig: str = ""

    def __post_init__(self) -> None:
        object.__setattr__(self, "position", _coerce_position(self.position, "position"))
        object.__setattr__(self, "read_count", _coerce_read_count(self.read_count))
        object.__setattr__(self, "contig", _coerce_contig(self.contig))


@dataclass(frozen=True)
class GeneAssignment:
    """The CDS to which a candidate insertion site belongs."""

    gene_id: str
    strand: str

    def __post_init__(self) -> None:
        object.__setattr__(self, "gene_id", _coerce_gene_id(self.gene_id))
        object.__setattr__(self, "strand", _coerce_strand(self.strand))


@dataclass(frozen=True)
class AnnotatedInsertionSite:
    """An insertion site annotated with zero, one, or several overlapping CDSs."""

    position: int
    read_count: float
    contig: str = ""
    genes: tuple[GeneAssignment, ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(self, "position", _coerce_position(self.position, "position"))
        object.__setattr__(self, "read_count", _coerce_read_count(self.read_count))
        object.__setattr__(self, "contig", _coerce_contig(self.contig))
        assignments = tuple(self.genes)
        if not all(isinstance(assignment, GeneAssignment) for assignment in assignments):
            raise TypeError("genes must contain GeneAssignment objects")
        object.__setattr__(self, "genes", assignments)


@dataclass(frozen=True)
class GeneSiteRow:
    """One gene-specific view of a site passed to a custom classifier.

    An overlapping site appears once for each overlapping gene, which is the
    same flattening performed by ``separate_rows(gene, sep = ",")`` in the R
    workflow.
    """

    contig: str
    gene_id: str
    strand: str
    position: int
    read_count: float


@dataclass(frozen=True)
class GapEvidence:
    """How surprising this gene's longest run of empty sites is.

    Saturation cannot see contiguity, so this is independent evidence rather
    than a restatement of the saturation-based call.  See
    :mod:`phage_tnseq_viz.gaps` for the model.

    ``qvalue`` is the Benjamini-Hochberg adjusted ``pvalue``; ``significant``
    records ``qvalue < GAP_ALPHA``.

    ``domain_candidate`` marks the disagreement worth looking at by hand: the run
    of dead sites is significant, yet the built-in rules did not call the gene
    essential.  A gene where one domain tolerates insertion and another does not
    produces exactly this pattern, and the saturation-based rules cannot see it
    because the gene's overall saturation looks ordinary.  It is a prompt to
    inspect the gene, not a call in its own right.
    """

    max_gap: int
    expected_max_gap: float
    pvalue: float
    qvalue: float
    significant: bool
    domain_candidate: bool


@dataclass(frozen=True)
class GeneClassification:
    """Per-gene statistics and the resulting essentiality label.

    ``initial_call`` records the saturation-only screen.  ``final_call`` is the
    R refinement/fitness-adjusted result, or a custom result when a custom
    classifier is supplied.  ``None`` faithfully represents R's ``NA`` for a
    gene with five or fewer candidate sites.

    ``read_count_median_threshold`` is the pooled positive-count median of *this
    gene's contig* used by the two fitness refinements; it is ``None`` when the
    contig had no positive sites.

    ``binomial_min_sites`` is the site count above which an empty gene on this
    contig is called ``Essential (binomial)``; it is ``None`` when the supplement
    is disabled or the contig's saturation admits no threshold.

    ``gap`` carries the independent run-length evidence, or ``None`` when the
    gap analysis is disabled.  It never affects ``final_call``.
    """

    contig: str
    gene_id: str
    strand: str
    total_sites: int
    hits: int
    saturation: float
    initial_call: str | None
    final_call: str | None
    read_count_median_threshold: float | None = None
    binomial_min_sites: float | None = None
    gap: GapEvidence | None = None


@dataclass(frozen=True)
class ClassificationResult:
    """The calls plus the positive-count median(s) used by the R method.

    ``read_thresholds`` maps each contig to its own pooled positive-count median.
    Thresholds are kept per contig because the fitness refinements are a
    within-genome comparison, and phage genome sequencing depths are not
    comparable across contigs.  ``read_threshold`` is retained as a convenience
    for the common single-contig dataset (the sole contig's threshold), and is
    ``None`` whenever more than one contig is present.

    ``library_saturation`` and ``binomial_min_sites`` record, per contig, the
    genic saturation and the site count derived from it by
    :func:`binomial_min_sites`.  Both are empty when the binomial supplement is
    disabled.  Saturation is worth reporting on its own: it is the number that
    decides how much the supplement can recover, and a very low value is the
    usual sign of an under-sequenced library.
    """

    calls: tuple[GeneClassification, ...]
    read_threshold: float | None
    read_thresholds: dict[str, float | None] = field(default_factory=dict)
    library_saturation: dict[str, float] = field(default_factory=dict)
    binomial_min_sites: dict[str, float | None] = field(default_factory=dict)


@dataclass(frozen=True)
class ConsensusCall:
    """A consensus label over subsampling replicates."""

    contig: str
    gene_id: str
    consensus_call: str
    votes: int
    n_replicates: int
    replicate_calls: tuple[str | None, ...]


GeneClassifier = Callable[[str, tuple[GeneSiteRow, ...]], Optional[str]]


@dataclass(frozen=True)
class _PendingCall:
    """One gene's result before gap q-values are known across the whole set."""

    call: GeneClassification
    called_essential: bool
    max_gap: int
    expected_gap: float
    pvalue: float


def annotate_sites_with_genes(
    sites: Iterable[InsertionSite],
    genes: Iterable[object],
    *,
    contig: str | None = None,
    gene_id_getter: Callable[[object], str | None] | None = None,
) -> list[AnnotatedInsertionSite]:
    """Annotate candidate sites using 1-based inclusive CDS coordinates.

    ``genes`` may be the project's :class:`~phage_tnseq_viz.genome.Gene`
    objects or mapping/object records with ``start``, ``end``, ``strand``, and
    optionally ``locus``.  The default identifier is ``gene.locus`` (currently
    populated from ``locus_tag`` before ``ID``/``label``), falling back to the
    deterministic ``contig:start-end:strand`` requested by the plotting/CSV
    layer.  Pass the record accession through ``contig=...`` whenever sites do
    not already carry it.

    Sites outside CDSs are retained with an empty ``genes`` tuple.  A site that
    overlaps several CDSs receives every assignment, as in the source R code.
    """

    input_sites = [_coerce_site(site) for site in sites]
    annotation_contig = _resolve_annotation_contig(input_sites, contig)

    prepared_genes: list[tuple[int, int, GeneAssignment]] = []
    for gene in genes:
        start = _coerce_position(_value(gene, "start"), "gene start")
        end = _coerce_position(_value(gene, "end"), "gene end")
        if end < start:
            raise ValueError(f"gene end ({end}) is before start ({start})")
        raw_strand = _value(gene, "strand")
        strand = _coerce_strand(raw_strand)
        supplied_id = gene_id_getter(gene) if gene_id_getter is not None else _optional_value(gene, "locus")
        if supplied_id is None or not str(supplied_id).strip():
            if not annotation_contig:
                raise ValueError(
                    "contig is required to create fallback gene IDs for genes without a locus tag"
                )
            # Keep this exactly compatible with genome.gene_identifier(), whose
            # ``Gene.strand`` is stored as +1/-1 rather than +/- symbols.
            supplied_id = f"{annotation_contig}:{start}-{end}:{raw_strand}"
        prepared_genes.append((start, end, GeneAssignment(str(supplied_id), strand)))

    annotated: list[AnnotatedInsertionSite] = []
    for site in input_sites:
        if annotation_contig and site.contig and site.contig != annotation_contig:
            raise ValueError(
                f"site contig {site.contig!r} does not match annotation contig {annotation_contig!r}"
            )
        site_contig = site.contig or annotation_contig
        assignments = tuple(
            assignment
            for start, end, assignment in prepared_genes
            if start <= site.position <= end
        )
        annotated.append(
            AnnotatedInsertionSite(
                position=site.position,
                read_count=site.read_count,
                contig=site_contig,
                genes=assignments,
            )
        )
    return annotated


def classify_genes(
    sites: Iterable[AnnotatedInsertionSite],
    *,
    classifier: GeneClassifier | None = None,
    binomial_short_genes: bool = True,
    gap_analysis: bool = True,
) -> ClassificationResult:
    """Classify genes using the audited R workflow or a supplied classifier.

    The read threshold is the R default ``quantile(positive_counts, 0.50)`` after
    sites are flattened across overlapping genes, computed *per contig*: the R
    script assumes a single genome, and pooling a shallow contig's counts with a
    deeper one's would shift the shared median and silently flip fitness calls.
    When ``classifier`` is provided, R statistics and ``initial_call`` are still
    calculated, but its return value replaces every ``final_call``.

    ``binomial_short_genes`` adds the :func:`binomial_min_sites` supplement for
    genes the R rules leave uncalled.  It only ever fills a ``None`` result, so
    it cannot change a call the R rules made; set it to ``False`` to reproduce
    the R script's output exactly.  Saturation is measured per contig, for the
    same reason the read thresholds are.

    ``gap_analysis`` attaches :class:`GapEvidence` to every gene.  It is purely
    additional evidence and never changes ``final_call``, so it is safe to leave
    on; the two methods answer different questions and are most useful where
    they disagree.  Its p-values are corrected together across every gene in the
    result, contigs included, because the multiple-testing family is the set of
    tests actually performed rather than a per-genome comparison.
    """

    groups = _group_sites(sites)
    positive_by_contig: dict[str, list[float]] = {}
    site_totals: dict[str, int] = {}
    hit_totals: dict[str, int] = {}
    for (contig, _gene_id), rows in groups.items():
        site_totals[contig] = site_totals.get(contig, 0) + len(rows)
        for row in rows:
            if row.read_count > 0:
                hit_totals[contig] = hit_totals.get(contig, 0) + 1
                positive_by_contig.setdefault(contig, []).append(row.read_count)
    read_thresholds: dict[str, float | None] = {
        contig: _r_quantile(counts, 0.50)
        for contig, counts in positive_by_contig.items()
    }

    # Saturation is summed over genes, matching how the reference implementation
    # of this supplement derives its library-wide insertion frequency.  Genes are
    # therefore the denominator, not the whole contig; a site shared by two
    # overlapping CDSs counts once per gene, as it does everywhere else here.
    library_saturation: dict[str, float] = {
        contig: (hit_totals.get(contig, 0) / total if total else 0.0)
        for contig, total in site_totals.items()
    }
    binomial_thresholds: dict[str, float | None] = {}
    if binomial_short_genes:
        binomial_thresholds = {
            contig: binomial_min_sites(saturation)
            for contig, saturation in library_saturation.items()
        }

    # Built in two passes: the gap q-values need every p-value in hand.
    pending: list[_PendingCall] = []
    for (contig, gene_id), rows in groups.items():
        read_threshold = read_thresholds.get(contig)
        strand = rows[0].strand  # R uses unique(gene_data$strand)[1].
        total_sites = len(rows)
        hits = sum(row.read_count > 0 for row in rows)
        saturation = hits / total_sites
        initial_call, final_call = _r_classify_gene(rows, read_threshold)

        threshold_sites = binomial_thresholds.get(contig)
        if (
            binomial_short_genes
            and final_call is None
            and hits == 0
            and threshold_sites is not None
            and total_sites > threshold_sites
        ):
            final_call = ESSENTIAL_BINOMIAL

        # Recorded before any custom classifier runs, so the gap comparison keeps
        # its meaning when a caller replaces the labels with its own vocabulary.
        called_essential = final_call in _ESSENTIAL_CALLS

        max_gap = 0
        pvalue = 1.0
        expected_gap = 0.0
        if gap_analysis:
            saturation_here = library_saturation.get(contig, 0.0)
            ordered = sorted(rows, key=lambda row: row.position)
            max_gap = gaps.longest_zero_run(row.read_count for row in ordered)
            pvalue = gaps.gap_pvalue(total_sites, max_gap, saturation_here)
            expected_gap = _expected_gap(total_sites, saturation_here)

        if classifier is not None:
            try:
                custom_call = classifier(gene_id, tuple(rows))
            except Exception as exc:  # noqa: BLE001 - add gene context to user code errors
                raise RuntimeError(
                    f"custom classifier failed for {contig}:{gene_id}"
                ) from exc
            final_call = _validate_custom_call(custom_call)

        pending.append(
            _PendingCall(
                call=GeneClassification(
                    contig=contig,
                    gene_id=gene_id,
                    strand=strand,
                    total_sites=total_sites,
                    hits=hits,
                    saturation=saturation,
                    initial_call=initial_call,
                    final_call=final_call,
                    read_count_median_threshold=read_threshold,
                    binomial_min_sites=threshold_sites,
                ),
                called_essential=called_essential,
                max_gap=max_gap,
                expected_gap=expected_gap,
                pvalue=pvalue,
            )
        )

    if gap_analysis:
        qvalues = gaps.benjamini_hochberg([entry.pvalue for entry in pending])
        calls = [
            replace(
                entry.call,
                gap=GapEvidence(
                    max_gap=entry.max_gap,
                    expected_max_gap=entry.expected_gap,
                    pvalue=entry.pvalue,
                    qvalue=qvalue,
                    significant=qvalue < GAP_ALPHA,
                    domain_candidate=qvalue < GAP_ALPHA and not entry.called_essential,
                ),
            )
            for entry, qvalue in zip(pending, qvalues)
        ]
    else:
        calls = [entry.call for entry in pending]

    contigs = {contig for contig, _gene_id in groups}
    single_threshold = read_thresholds.get(next(iter(contigs))) if len(contigs) == 1 else None
    return ClassificationResult(
        calls=tuple(calls),
        read_threshold=single_threshold,
        read_thresholds=read_thresholds,
        library_saturation=library_saturation,
        binomial_min_sites=binomial_thresholds,
    )


def _expected_gap(total_sites: int, saturation: float) -> float:
    """Expected longest empty run, or 0.0 when the Bernoulli model does not apply.

    A contig where every candidate site was hit, or none was, has no meaningful
    expected run; both degenerate cases also give a p-value of 1.0, so reporting
    zero here keeps the two consistent.
    """

    non_insertion = 1.0 - saturation
    if total_sites <= 0 or not 0.0 < non_insertion < 1.0:
        return 0.0
    return gaps.expected_longest_run(total_sites, non_insertion)


def load_classifier(path: str | Path) -> GeneClassifier:
    """Load ``classify_gene(gene_id, site_rows)`` from a trusted local ``.py`` file.

    The file is intentionally loaded as ordinary Python, so it has the same
    privileges as the process running this application.
    """

    source = Path(path).expanduser().resolve()
    if not source.is_file():
        raise FileNotFoundError(f"custom classifier file not found: {source}")
    if source.suffix.lower() != ".py":
        raise ValueError(f"custom classifier must be a .py file: {source}")

    digest = hashlib.sha256(str(source).encode("utf-8")).hexdigest()[:16]
    module_name = f"_phage_tnseq_custom_classifier_{digest}"
    spec = importlib.util.spec_from_file_location(module_name, source)
    if spec is None or spec.loader is None:  # pragma: no cover - defensive importlib guard
        raise ImportError(f"could not create an import specification for {source}")
    module = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(module)
    except Exception as exc:  # noqa: BLE001 - show the source file in the error
        raise RuntimeError(f"could not load custom classifier {source}") from exc

    classifier = getattr(module, "classify_gene", None)
    if not callable(classifier):
        raise ValueError(
            f"custom classifier {source} must define callable "
            "classify_gene(gene_id, site_rows)"
        )
    return classifier


def apply_classifier(
    sites: Iterable[AnnotatedInsertionSite],
    classifier: GeneClassifier,
) -> ClassificationResult:
    """Apply a loaded/custom classifier while retaining standard site statistics."""

    return classify_genes(sites, classifier=classifier)


def consensus_calls(
    replicates: Sequence[Sequence[GeneClassification] | ClassificationResult],
    *,
    min_votes: int = 3,
) -> tuple[ConsensusCall, ...]:
    """Return deterministic consensus calls across subsampling replicates.

    The supplied R script assumes five replicates and accepts a label seen at
    least three times.  Those are the defaults here.  Missing/``None`` calls do
    not count as a label, and a tie is reported as ``Ambiguous`` rather than
    inheriting input-order behaviour from R's ``slice_head``.
    """

    if min_votes < 1:
        raise ValueError("min_votes must be at least 1")

    normalised: list[dict[tuple[str, str], str | None]] = []
    all_keys: set[tuple[str, str]] = set()
    for replicate_index, replicate in enumerate(replicates, start=1):
        source_calls = replicate.calls if isinstance(replicate, ClassificationResult) else replicate
        by_gene: dict[tuple[str, str], str | None] = {}
        for call in source_calls:
            if not isinstance(call, GeneClassification):
                raise TypeError(
                    "replicates must contain GeneClassification sequences or ClassificationResult objects"
                )
            key = (call.contig, call.gene_id)
            if key in by_gene:
                raise ValueError(
                    f"replicate {replicate_index} contains duplicate gene {call.contig}:{call.gene_id}"
                )
            by_gene[key] = call.final_call
            all_keys.add(key)
        normalised.append(by_gene)

    results: list[ConsensusCall] = []
    for contig, gene_id in sorted(all_keys):
        per_replicate = tuple(replicate.get((contig, gene_id)) for replicate in normalised)
        votes = Counter(call for call in per_replicate if call is not None)
        if not votes:
            consensus = AMBIGUOUS
            support = 0
        else:
            support = max(votes.values())
            leaders = [label for label, count in votes.items() if count == support]
            consensus = leaders[0] if len(leaders) == 1 and support >= min_votes else AMBIGUOUS
        results.append(
            ConsensusCall(
                contig=contig,
                gene_id=gene_id,
                consensus_call=consensus,
                votes=support,
                n_replicates=len(normalised),
                replicate_calls=per_replicate,
            )
        )
    return tuple(results)


def _group_sites(
    sites: Iterable[AnnotatedInsertionSite],
) -> dict[tuple[str, str], list[GeneSiteRow]]:
    """Flatten overlapping gene annotations exactly as the R input does."""

    groups: dict[tuple[str, str], list[GeneSiteRow]] = {}
    for site in sites:
        annotated = _coerce_annotated_site(site)
        for assignment in annotated.genes:
            key = (annotated.contig, assignment.gene_id)
            groups.setdefault(key, []).append(
                GeneSiteRow(
                    contig=annotated.contig,
                    gene_id=assignment.gene_id,
                    strand=assignment.strand,
                    position=annotated.position,
                    read_count=annotated.read_count,
                )
            )
    return groups


def _r_classify_gene(
    rows: Sequence[GeneSiteRow], read_threshold: float | None
) -> tuple[str | None, str | None]:
    """Port one gene's classification from the supplied R script."""

    total_sites = len(rows)
    hits = sum(row.read_count > 0 for row in rows)
    saturation = hits / total_sites

    if total_sites <= 5:
        return None, None
    if saturation < 0.2:
        return ESSENTIAL, ESSENTIAL
    if saturation > 0.8:
        return NON_ESSENTIAL, NON_ESSENTIAL

    # ``Further analysis`` is R's temporary label; every qualifying gene then
    # receives either Essential or Intermediate from classify_gene().
    initial_call = FURTHER_ANALYSIS
    strand = rows[0].strand
    hit_positions = [row.position for row in rows if row.read_count > 0]
    positions = [row.position for row in rows]
    position_threshold = _r_quantile(positions, 0.80 if strand == "+" else 0.20)
    median_hit_position = _r_quantile(hit_positions, 0.50)

    if (
        strand == "+" and median_hit_position is not None and position_threshold is not None
        and median_hit_position >= position_threshold
    ):
        final_call = ESSENTIAL
    elif (
        strand == "-" and median_hit_position is not None and position_threshold is not None
        and median_hit_position <= position_threshold
    ):
        final_call = ESSENTIAL
    else:
        oriented = sorted(rows, key=lambda row: row.position, reverse=strand == "-")
        midpoint = total_sites // 2
        first_half = oriented[:midpoint]
        second_half = oriented[midpoint:]
        zeros_first = sum(row.read_count == 0 for row in first_half) / len(first_half)
        zeros_second = sum(row.read_count == 0 for row in second_half) / len(second_half)
        if (zeros_first > 0.7 and zeros_second < 0.3) or (
            zeros_second > 0.7 and zeros_first < 0.3
        ):
            final_call = ESSENTIAL
        else:
            final_call = INTERMEDIATE

    # The two post-refinement adjustments mirror the R script's strict count
    # comparisons and inclusive saturation windows.
    positive_counts = [row.read_count for row in rows if row.read_count > 0]
    if read_threshold is not None and 0.2 <= saturation <= 0.4:
        prop_weak = sum(count < read_threshold for count in positive_counts) / len(positive_counts)
        if prop_weak > 0.5 and final_call != ESSENTIAL:
            final_call = STRONG_FITNESS_DEFECT
    if read_threshold is not None and 0.7 <= saturation <= 0.8:
        prop_strong = sum(count > read_threshold for count in positive_counts) / len(positive_counts)
        if prop_strong > 0.5 and final_call != ESSENTIAL:
            final_call = REDUCED_FITNESS

    return initial_call, final_call


def _r_quantile(values: Sequence[float], probability: float) -> float | None:
    """R's default type-7 quantile, used by the supplied script."""

    if not values:
        return None
    ordered = sorted(float(value) for value in values)
    h = (len(ordered) - 1) * probability
    low = math.floor(h)
    high = math.ceil(h)
    if low == high:
        return ordered[low]
    return ordered[low] + (h - low) * (ordered[high] - ordered[low])


def _coerce_site(site: object) -> InsertionSite:
    if isinstance(site, InsertionSite):
        return site
    return InsertionSite(
        position=_value(site, "position"),
        read_count=_value(site, "read_count"),
        contig=_optional_value(site, "contig") or "",
    )


def _coerce_annotated_site(site: object) -> AnnotatedInsertionSite:
    if isinstance(site, AnnotatedInsertionSite):
        return site
    return AnnotatedInsertionSite(
        position=_value(site, "position"),
        read_count=_value(site, "read_count"),
        contig=_optional_value(site, "contig") or "",
        genes=tuple(_value(site, "genes")),
    )


def _resolve_annotation_contig(sites: Sequence[InsertionSite], contig: str | None) -> str:
    if contig is not None:
        return _coerce_contig(contig)
    present = {site.contig for site in sites if site.contig}
    if len(present) > 1:
        raise ValueError(
            "annotate one contig at a time, or pass a single matching contig explicitly"
        )
    return next(iter(present), "")


def _value(record: object, name: str) -> Any:
    value = _optional_value(record, name)
    if value is None:
        raise ValueError(f"record is missing required field {name!r}")
    return value


def _optional_value(record: object, name: str) -> Any | None:
    if isinstance(record, Mapping):
        return record.get(name)
    return getattr(record, name, None)


def _coerce_position(value: object, label: str) -> int:
    if isinstance(value, bool):
        raise ValueError(f"{label} must be a positive integer, not bool")
    try:
        numeric = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} must be a positive integer: {value!r}") from exc
    if not math.isfinite(numeric) or not numeric.is_integer() or numeric < 1:
        raise ValueError(f"{label} must be a positive integer: {value!r}")
    return int(numeric)


def _coerce_read_count(value: object) -> float:
    if isinstance(value, bool):
        raise ValueError("read_count must be a non-negative number, not bool")
    try:
        count = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"read_count must be a non-negative number: {value!r}") from exc
    if not math.isfinite(count) or count < 0:
        raise ValueError(f"read_count must be a non-negative number: {value!r}")
    return count


def _coerce_contig(value: object) -> str:
    if value is None:
        return ""
    return str(value).strip()


def _coerce_gene_id(value: object) -> str:
    if value is None or not str(value).strip():
        raise ValueError("gene_id must be a non-empty string")
    return str(value).strip()


def _coerce_strand(value: object) -> str:
    if value in (1, "+1"):
        return "+"
    if value in (-1, "-1"):
        return "-"
    text = str(value).strip()
    if text == "+":
        return "+"
    if text == "-":
        return "-"
    raise ValueError(f"strand must be '+', '-', +1, or -1: {value!r}")


def _validate_custom_call(value: object) -> str | None:
    if value is None:
        return None
    if not isinstance(value, str) or not value.strip():
        raise ValueError("custom classify_gene() must return a non-empty str or None")
    return value.strip()


__all__ = [
    "AMBIGUOUS",
    "BINOMIAL_ALPHA",
    "ESSENTIAL",
    "ESSENTIAL_BINOMIAL",
    "FURTHER_ANALYSIS",
    "GAP_ALPHA",
    "INTERMEDIATE",
    "NON_ESSENTIAL",
    "REDUCED_FITNESS",
    "STRONG_FITNESS_DEFECT",
    "AnnotatedInsertionSite",
    "ClassificationResult",
    "ConsensusCall",
    "GapEvidence",
    "GeneAssignment",
    "GeneClassification",
    "GeneClassifier",
    "GeneSiteRow",
    "InsertionSite",
    "annotate_sites_with_genes",
    "apply_classifier",
    "binomial_min_sites",
    "classify_genes",
    "consensus_calls",
    "load_classifier",
]
