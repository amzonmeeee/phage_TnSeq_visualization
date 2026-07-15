"""Generate small synthetic phage GenBank files and a final Tn-Seq table.

Creates:
  * example_annotated.gbk   – with phold-style 'function' qualifiers
  * example_bare.gbk        – sequence only, no CDS annotation
  * example_final_sites.csv – complete synthetic Himar1/TA final dataset

The CSV deliberately contains every TA candidate site, including zeros.  It is
therefore useful for testing the fast ``plot --final-dataset`` path without
installing FastQC/fastp/TRANSIT/BWA/SeqKit or fabricating a protocol-specific
raw FASTQ library.
"""

from __future__ import annotations

import random
import csv
import argparse
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

random.seed(42)

CATEGORIES = [
    "head and packaging", "tail", "connector", "lysis",
    "dna, rna and nucleotide metabolism", "integration and excision",
    "transcription regulation", "moron, auxiliary metabolic gene and host takeover",
    "other", "unknown function",
]


def random_seq(n: int) -> str:
    return "".join(random.choice("ACGT") for _ in range(n))


def build(annotated: bool) -> SeqRecord:
    seq = random_seq(20000)
    rec = SeqRecord(Seq(seq), id="NC_TEST01", name="phiTEST",
                    description="Synthetic test phage phiTEST")
    rec.annotations["molecule_type"] = "DNA"
    rec.annotations["topology"] = "linear"

    if annotated:
        pos = 200
        i = 0
        while pos < 19000:
            length = random.choice([300, 450, 600, 900, 1200])
            strand = random.choice([1, -1])
            cat = CATEGORIES[i % len(CATEGORIES)]
            feat = SeqFeature(
                FeatureLocation(pos, pos + length, strand=strand),
                type="CDS",
                qualifiers={
                    "function": [cat],
                    "product": [f"{cat} protein"],
                    "locus_tag": [f"phiTEST_{i+1:03d}"],
                },
            )
            rec.features.append(feat)
            pos += length + random.randint(80, 300)
            i += 1
        # a couple of tRNAs
        for tpos in (5000, 12000):
            rec.features.append(SeqFeature(
                FeatureLocation(tpos, tpos + 75, strand=1),
                type="tRNA", qualifiers={"product": ["tRNA-Leu"]},
            ))
    return rec


def write_final_sites(record: SeqRecord, path: Path) -> None:
    """Create reproducible TA-site counts with several fitness behaviours."""
    rng = random.Random(3900)
    sequence = str(record.seq).upper()
    candidates = [position + 1 for position in range(len(sequence) - 1)
                  if sequence[position:position + 2] == "TA"]
    cds = [feature for feature in record.features if feature.type == "CDS"]

    def gene_index(position: int) -> int | None:
        for index, feature in enumerate(cds):
            if int(feature.location.start) < position <= int(feature.location.end):
                return index
        return None

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("contig", "position", "read_count"))
        writer.writeheader()
        for position in candidates:
            index = gene_index(position)
            # Cycle CDSs through saturated/non-essential, depleted/essential,
            # and intermediate patterns. Intergenic positions are modestly hit.
            pattern = index % 3 if index is not None else 2
            hit_probability = (0.92, 0.08, 0.48)[pattern]
            if rng.random() < hit_probability:
                read_count = rng.randint(8, 35) if pattern == 0 else rng.randint(1, 16)
            else:
                read_count = 0
            writer.writerow({"contig": record.id, "position": position, "read_count": read_count})


def main(argv: list[str] | None = None) -> None:
    from Bio import SeqIO
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).parent)
    args = parser.parse_args(argv)
    out_dir = args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    annotated = build(True)
    bare = build(False)
    SeqIO.write(annotated, out_dir / "example_annotated.gbk", "genbank")
    SeqIO.write(bare, out_dir / "example_bare.gbk", "genbank")
    write_final_sites(annotated, out_dir / "example_final_sites.csv")
    print("Wrote example_annotated.gbk, example_bare.gbk, and example_final_sites.csv")


if __name__ == "__main__":
    main()
