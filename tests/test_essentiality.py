"""Focused tests for the dependency-free R essentiality-classifier port."""

from __future__ import annotations

from dataclasses import dataclass

import pytest

from phage_tnseq_viz.essentiality import (
    AMBIGUOUS,
    ESSENTIAL,
    FURTHER_ANALYSIS,
    INTERMEDIATE,
    NON_ESSENTIAL,
    REDUCED_FITNESS,
    STRONG_FITNESS_DEFECT,
    AnnotatedInsertionSite,
    GeneAssignment,
    GeneClassification,
    InsertionSite,
    annotate_sites_with_genes,
    apply_classifier,
    classify_genes,
    consensus_calls,
    load_classifier,
)
from phage_tnseq_viz.genome import gene_identifier


@dataclass
class _Gene:
    start: int
    end: int
    strand: int
    locus: str | None = None


def _annotated_rows(
    gene_id: str,
    strand: str,
    counts: list[float],
    *,
    start: int = 1,
    contig: str = "ctg",
) -> list[AnnotatedInsertionSite]:
    """Make consecutive, fully annotated candidate-site rows for one CDS."""

    assignment = GeneAssignment(gene_id, strand)
    return [
        AnnotatedInsertionSite(
            position=start + offset,
            read_count=count,
            contig=contig,
            genes=(assignment,),
        )
        for offset, count in enumerate(counts)
    ]


def _calls_by_gene(rows: list[AnnotatedInsertionSite]):
    return {call.gene_id: call for call in classify_genes(rows).calls}


def _call(gene_id: str, final_call: str | None, *, contig: str = "ctg") -> GeneClassification:
    return GeneClassification(
        contig=contig,
        gene_id=gene_id,
        strand="+",
        total_sites=10,
        hits=5,
        saturation=0.5,
        initial_call=FURTHER_ANALYSIS,
        final_call=final_call,
    )


def test_annotation_keeps_intergenic_sites_and_expands_overlaps() -> None:
    sites = [
        InsertionSite(position=10, read_count=4, contig="NC_TEST"),
        InsertionSite(position=20, read_count=0, contig="NC_TEST"),
        InsertionSite(position=30, read_count=2, contig="NC_TEST"),
    ]
    genes = [
        _Gene(start=5, end=20, strand=1, locus="known_locus"),
        _Gene(start=20, end=25, strand=-1),
    ]

    annotated = annotate_sites_with_genes(sites, genes, contig="NC_TEST")

    assert [assignment.gene_id for assignment in annotated[0].genes] == ["known_locus"]
    assert [(assignment.gene_id, assignment.strand) for assignment in annotated[1].genes] == [
        ("known_locus", "+"),
        ("NC_TEST:20-25:-1", "-"),
    ]
    assert annotated[1].genes[1].gene_id == gene_identifier(genes[1], contig="NC_TEST")
    assert annotated[2].genes == ()


def test_screening_thresholds_and_small_genes_match_r() -> None:
    rows = (
        _annotated_rows("low", "+", [0] * 6, start=1)
        + _annotated_rows("high", "+", [5] * 6, start=20)
        + _annotated_rows("small", "+", [0] * 5, start=40)
    )

    calls = _calls_by_gene(rows)

    assert calls["low"].saturation == 0
    assert calls["low"].initial_call == ESSENTIAL
    assert calls["low"].final_call == ESSENTIAL
    assert calls["high"].saturation == 1
    assert calls["high"].initial_call == NON_ESSENTIAL
    assert calls["high"].final_call == NON_ESSENTIAL
    assert calls["small"].initial_call is None
    assert calls["small"].final_call is None


def test_three_prime_bias_is_strand_aware() -> None:
    # At saturation exactly 0.2 these remain in R's "Further analysis" band.
    # + strand hits cluster at the high-coordinate 3' end; - strand hits cluster
    # at the low-coordinate 3' end.
    rows = (
        _annotated_rows("plus", "+", [0] * 8 + [1, 1], start=1)
        + _annotated_rows("minus", "-", [1, 1] + [0] * 8, start=101)
    )

    calls = _calls_by_gene(rows)

    assert calls["plus"].initial_call == FURTHER_ANALYSIS
    assert calls["plus"].final_call == ESSENTIAL
    assert calls["minus"].initial_call == FURTHER_ANALYSIS
    assert calls["minus"].final_call == ESSENTIAL


def test_domain_essentiality_and_intermediate_classification() -> None:
    # The domain case has 80% zero sites in the first oriented half and 20% in
    # the second, but its hit median is below the 80th percentile so it reaches
    # the second R refinement rule.
    domain_counts = [0, 0, 0, 1, 0, 1, 1, 1, 1, 0]
    intermediate_counts = [0, 0, 0, 1, 1, 1, 1, 1, 1, 0]
    rows = (
        _annotated_rows("domain", "+", domain_counts, start=1)
        + _annotated_rows("intermediate", "+", intermediate_counts, start=101)
    )

    calls = _calls_by_gene(rows)

    assert calls["domain"].final_call == ESSENTIAL
    assert calls["intermediate"].final_call == INTERMEDIATE


def test_borderline_fitness_adjustments_use_pooled_positive_median() -> None:
    # The all-hit reference sites make the pooled positive median 100.  This
    # leaves low_saturation as an intermediate call before its weak-read
    # adjustment, and high_saturation as intermediate before its strong-read
    # adjustment.
    rows = (
        _annotated_rows("low_saturation", "+", [0, 1, 0, 1, 0, 1, 0, 0, 0, 0], start=1)
        + _annotated_rows("reference", "+", [100] * 10, start=101)
        + _annotated_rows("high_saturation", "+", [200] * 7 + [0] * 3, start=201)
        # This one is three-prime Essential at 0.2 saturation and must *not*
        # be replaced by Strong fitness defect despite weak positive counts.
        + _annotated_rows("protected_essential", "+", [0] * 8 + [1, 1], start=301)
    )

    result = classify_genes(rows)
    calls = {call.gene_id: call for call in result.calls}

    assert result.read_threshold == pytest.approx(100)
    assert calls["low_saturation"].final_call == STRONG_FITNESS_DEFECT
    assert calls["high_saturation"].final_call == REDUCED_FITNESS
    assert calls["protected_essential"].final_call == ESSENTIAL


def test_read_threshold_is_computed_per_contig() -> None:
    # The two fitness refinements are a within-genome comparison, so each contig
    # must use its own pooled positive-count median.  Contig A (deep) alone lands
    # on Intermediate at saturation 0.75; pooling a shallow contig B used to drag
    # the shared median down far enough that A's counts all read as "strong" and
    # it flipped to Reduced fitness.  With per-contig thresholds A is stable.
    deep = _annotated_rows("geneA", "+", [1000, 900, 1100, 950, 1050, 1200, 0, 0], contig="A")
    shallow = _annotated_rows("geneB", "+", [5, 6, 4, 7, 5, 6, 0, 0], contig="B")

    alone = {call.gene_id: call for call in classify_genes(deep).calls}
    pooled_result = classify_genes(deep + shallow)
    pooled = {call.gene_id: call for call in pooled_result.calls}

    assert alone["geneA"].final_call == INTERMEDIATE
    assert pooled["geneA"].final_call == INTERMEDIATE
    # Each gene reports its own contig's median, unaffected by the other contig.
    assert pooled["geneA"].read_count_median_threshold == alone["geneA"].read_count_median_threshold
    assert pooled_result.read_thresholds["A"] != pooled_result.read_thresholds["B"]
    # The single-value convenience is None once more than one contig is present.
    assert pooled_result.read_threshold is None
    assert classify_genes(deep).read_threshold == pooled_result.read_thresholds["A"]


def test_consensus_requires_three_votes_and_makes_ties_ambiguous() -> None:
    replicates = [
        [_call("stable", ESSENTIAL), _call("tied", ESSENTIAL), _call("missing", None)],
        [_call("stable", ESSENTIAL), _call("tied", NON_ESSENTIAL), _call("missing", ESSENTIAL)],
        [_call("stable", ESSENTIAL), _call("tied", ESSENTIAL), _call("missing", ESSENTIAL)],
        [_call("stable", NON_ESSENTIAL), _call("tied", NON_ESSENTIAL), _call("missing", ESSENTIAL)],
    ]

    calls = {call.gene_id: call for call in consensus_calls(replicates)}

    assert calls["stable"].consensus_call == ESSENTIAL
    assert calls["stable"].votes == 3
    assert calls["tied"].consensus_call == AMBIGUOUS
    assert calls["tied"].votes == 2
    assert calls["missing"].consensus_call == ESSENTIAL
    assert calls["missing"].replicate_calls[0] is None


def test_custom_classifier_loader_and_application(tmp_path) -> None:
    classifier_path = tmp_path / "my_classifier.py"
    classifier_path.write_text(
        "def classify_gene(gene_id, site_rows):\n"
        "    assert all(hasattr(row, 'position') for row in site_rows)\n"
        "    return 'Custom high' if sum(row.read_count for row in site_rows) >= 2 else None\n",
        encoding="utf-8",
    )
    classifier = load_classifier(classifier_path)
    rows = _annotated_rows("custom", "+", [0, 2, 1] + [0] * 7, start=1)

    result = apply_classifier(rows, classifier)

    assert result.calls[0].final_call == "Custom high"
    assert result.calls[0].initial_call == FURTHER_ANALYSIS


def test_custom_classifier_requires_the_documented_callable(tmp_path) -> None:
    bad_path = tmp_path / "bad_classifier.py"
    bad_path.write_text("def another_name():\n    return 'nope'\n", encoding="utf-8")

    with pytest.raises(ValueError, match="classify_gene"):
        load_classifier(bad_path)
