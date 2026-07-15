"""Template for a user-defined gene essentiality classifier.

Use with:

    phage-tnseq-viz plot genome.gbk --final-dataset final_sites.csv \
        --classifier examples/custom_classifier_template.py

``site_rows`` is a tuple of immutable ``GeneSiteRow`` objects.  Each object has
``contig``, ``gene_id``, ``strand``, ``position`` (1-based), and ``read_count``.
It includes zero-count candidate sites, and each overlap is supplied to every
overlapping CDS.  Return a short label string or ``None`` for no call.
"""

from __future__ import annotations


def classify_gene(gene_id, site_rows):
    """Return a custom call for one gene.

    Replace this deliberately simple example with your validated biological
    rule.  The label will be written verbatim to the gene CSV; unfamiliar labels
    are shown in neutral grey on the built-in map.
    """
    if len(site_rows) <= 5:
        return None

    saturation = sum(row.read_count > 0 for row in site_rows) / len(site_rows)
    if saturation < 0.20:
        return "Custom essential"
    if saturation > 0.80:
        return "Custom non-essential"
    return "Custom intermediate"
