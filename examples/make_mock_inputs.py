#!/usr/bin/env python3
"""Generate a matched mock final-dataset CSV and paired FASTQ files.

The reference GenBank sequence is the single source of truth.  Every Himar1/
mariner ``TA`` candidate is written to ``final-dataset.csv``, including sites
with zero reads.  Each non-zero CSV count is represented by exactly that many
read pairs in ``mock_reads1.fastq.gz`` and ``mock_reads2.fastq.gz``.

R1 models a transposon-junction read: it starts with a configurable terminal
primer and continues into the reference at a TA site.  R2 is taken from the
opposite end of the same genomic fragment.  These files are synthetic smoke-
test inputs, not a substitute for validating a real library's construct,
primer, error profile, or fragment-size distribution.

Example:

    python examples/make_mock_inputs.py /path/to/reference.gbk \
        --output-dir mock_inputs --read-pairs 5000 --seed 3900
"""

from __future__ import annotations

import argparse
from collections import Counter
from contextlib import contextmanager
import csv
from dataclasses import dataclass
import gzip
import hashlib
import io
import json
from pathlib import Path
import random
from typing import Iterator, Sequence, TextIO

from Bio import SeqIO
from Bio.Seq import reverse_complement
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord


DEFAULT_PRIMER = "ACTTATCAGCCAACCTGTTA"
GENE_HIT_PROBABILITIES = (0.92, 0.08, 0.48, 0.74, 0.28)
INTERGENIC_HIT_PROBABILITY = 0.55


@dataclass(frozen=True)
class CandidateSite:
    """One 1-based TA coordinate tied to its source record."""

    record_index: int
    contig: str
    position: int
    cds_index: int | None


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("reference", type=Path, help="Annotated GenBank reference.")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("mock_inputs"),
        help="Directory for the generated CSV, FASTQ files, and manifest.",
    )
    parser.add_argument(
        "--read-pairs",
        type=int,
        default=5000,
        help="Exact total number of paired reads to generate.",
    )
    parser.add_argument("--read-length", type=int, default=150)
    parser.add_argument("--insert-size", type=int, default=300)
    parser.add_argument(
        "--insert-jitter",
        type=int,
        default=30,
        help="Uniform plus/minus variation around --insert-size.",
    )
    parser.add_argument(
        "--primer",
        default=DEFAULT_PRIMER,
        help="Synthetic transposon terminal sequence placed at the start of R1.",
    )
    parser.add_argument("--seed", type=int, default=3900)
    return parser.parse_args(argv)


def load_records(path: Path) -> list[SeqRecord]:
    if not path.is_file():
        raise FileNotFoundError(f"GenBank reference not found: {path}")
    records = list(SeqIO.parse(str(path), "genbank"))
    if not records:
        raise ValueError(f"No GenBank records found in {path}")
    seen: set[str] = set()
    for record in records:
        record.id = record.id.strip()
        if not record.id:
            raise ValueError("Every GenBank record must have a non-blank ID")
        if record.id in seen:
            raise ValueError(f"Duplicate GenBank record ID: {record.id}")
        seen.add(record.id)
        sequence = str(record.seq).upper()
        if not sequence or set(sequence) - set("ACGTN"):
            raise ValueError(
                f"Record {record.id!r} must contain only A/C/G/T/N sequence"
            )
    return records


def enumerate_candidates(records: Sequence[SeqRecord]) -> list[CandidateSite]:
    candidates: list[CandidateSite] = []
    for record_index, record in enumerate(records):
        sequence = str(record.seq).upper()
        cds_features = [feature for feature in record.features if feature.type == "CDS"]
        for offset in range(len(sequence) - 1):
            if sequence[offset : offset + 2] != "TA":
                continue
            candidates.append(
                CandidateSite(
                    record_index=record_index,
                    contig=record.id,
                    position=offset + 1,
                    cds_index=_first_containing_cds(cds_features, offset),
                )
            )
    if not candidates:
        raise ValueError("The reference contains no TA candidate insertion sites")
    return candidates


def _first_containing_cds(
    cds_features: Sequence[SeqFeature], position_zero_based: int
) -> int | None:
    for index, feature in enumerate(cds_features):
        if position_zero_based in feature.location:
            return index
    return None


def _hit_probability(site: CandidateSite) -> float:
    if site.cds_index is None:
        return INTERGENIC_HIT_PROBABILITY
    return GENE_HIT_PROBABILITIES[site.cds_index % len(GENE_HIT_PROBABILITIES)]


def allocate_counts(
    candidates: Sequence[CandidateSite],
    records: Sequence[SeqRecord],
    *,
    read_pairs: int,
    minimum_fragment_length: int,
    rng: random.Random,
) -> Counter[CandidateSite]:
    """Choose positive sites and distribute exactly ``read_pairs`` among them."""
    counts: Counter[CandidateSite] = Counter()
    if read_pairs == 0:
        return counts

    eligible = [
        site
        for site in candidates
        if _maximum_available_fragment(site, records) >= minimum_fragment_length
    ]
    if not eligible:
        raise ValueError(
            "No TA site has enough flanking sequence for the requested read/insert size"
        )

    positive = [site for site in eligible if rng.random() < _hit_probability(site)]
    if not positive:
        positive = [rng.choice(eligible)]
    if len(positive) > read_pairs:
        positive = rng.sample(positive, read_pairs)

    counts.update(positive)
    remaining = read_pairs - len(positive)
    if remaining:
        # A log-normal weight creates ordinary depth variation and a few visibly
        # deeper sites without making the output depend on NumPy.
        weights = [rng.lognormvariate(0.0, 0.85) for _ in positive]
        counts.update(rng.choices(positive, weights=weights, k=remaining))
    return counts


def _maximum_available_fragment(
    site: CandidateSite, records: Sequence[SeqRecord]
) -> int:
    sequence_length = len(records[site.record_index].seq)
    ta_start = site.position - 1
    ta_end = ta_start + 2
    return max(sequence_length - ta_start, ta_end)


def write_final_dataset(
    path: Path,
    candidates: Sequence[CandidateSite],
    counts: Counter[CandidateSite],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=("contig", "position", "read_count")
        )
        writer.writeheader()
        for site in candidates:
            writer.writerow(
                {
                    "contig": site.contig,
                    "position": site.position,
                    "read_count": counts[site],
                }
            )


def write_fastq_pair(
    read1_path: Path,
    read2_path: Path,
    candidates: Sequence[CandidateSite],
    counts: Counter[CandidateSite],
    records: Sequence[SeqRecord],
    *,
    primer: str,
    read_length: int,
    insert_size: int,
    insert_jitter: int,
    rng: random.Random,
) -> int:
    genomic_r1_length = read_length - len(primer)
    minimum_fragment = insert_size - insert_jitter
    maximum_fragment = insert_size + insert_jitter
    read_number = 0
    with _open_deterministic_gzip_text(read1_path) as read1, \
            _open_deterministic_gzip_text(read2_path) as read2:
        for site in candidates:
            for _ in range(counts[site]):
                sequence = str(records[site.record_index].seq).upper()
                available = _maximum_available_fragment(site, records)
                fragment_length = rng.randint(
                    minimum_fragment, min(maximum_fragment, available)
                )
                read1_sequence, read2_sequence, strand = _make_pair_sequences(
                    sequence,
                    site.position,
                    fragment_length=fragment_length,
                    read_length=read_length,
                    genomic_r1_length=genomic_r1_length,
                    primer=primer,
                    rng=rng,
                )
                read_number += 1
                base_name = (
                    f"MOCK{read_number:08d}|contig={site.contig}"
                    f"|site={site.position}|strand={strand}"
                )
                _write_fastq_record(read1, f"{base_name}/1", read1_sequence)
                _write_fastq_record(read2, f"{base_name}/2", read2_sequence)
    return read_number


def _make_pair_sequences(
    sequence: str,
    position: int,
    *,
    fragment_length: int,
    read_length: int,
    genomic_r1_length: int,
    primer: str,
    rng: random.Random,
) -> tuple[str, str, str]:
    ta_start = position - 1
    ta_end = ta_start + 2
    orientations: list[str] = []
    if ta_start + fragment_length <= len(sequence):
        orientations.append("+")
    if ta_end - fragment_length >= 0:
        orientations.append("-")
    if not orientations:
        raise ValueError(
            f"TA site {position} lacks a {fragment_length}-bp genomic flank"
        )

    strand = rng.choice(orientations)
    if strand == "+":
        fragment_start = ta_start
        fragment_end = fragment_start + fragment_length
        junction_flank = sequence[
            fragment_start : fragment_start + genomic_r1_length
        ]
        opposite_end = reverse_complement(
            sequence[fragment_end - read_length : fragment_end]
        )
    else:
        fragment_end = ta_end
        fragment_start = fragment_end - fragment_length
        junction_flank = reverse_complement(
            sequence[fragment_end - genomic_r1_length : fragment_end]
        )
        opposite_end = sequence[fragment_start : fragment_start + read_length]

    read1_sequence = primer + str(junction_flank)
    read2_sequence = str(opposite_end)
    if len(read1_sequence) != read_length or len(read2_sequence) != read_length:
        raise AssertionError("Internal error: generated FASTQ read has the wrong length")
    if read1_sequence[len(primer) : len(primer) + 2] != "TA":
        raise AssertionError("Internal error: R1 genomic flank does not start at TA")
    return read1_sequence, read2_sequence, strand


def _write_fastq_record(handle: TextIO, name: str, sequence: str) -> None:
    quality = "I" * len(sequence)
    handle.write(f"@{name}\n{sequence}\n+\n{quality}\n")


@contextmanager
def _open_deterministic_gzip_text(path: Path) -> Iterator[TextIO]:
    """Write gzip without timestamps or path-dependent filename metadata."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as raw_handle:
        with gzip.GzipFile(
            filename="", mode="wb", fileobj=raw_handle, mtime=0
        ) as gzip_handle:
            with io.TextIOWrapper(
                gzip_handle, encoding="ascii", newline="\n"
            ) as text_handle:
                yield text_handle


def validate_args(args: argparse.Namespace) -> str:
    primer = args.primer.strip().upper()
    if set(primer) - set("ACGT") or not primer:
        raise ValueError("--primer must be a non-empty A/C/G/T sequence")
    if args.read_pairs < 0:
        raise ValueError("--read-pairs must be zero or greater")
    if args.read_length <= len(primer):
        raise ValueError("--read-length must be longer than the primer")
    if args.insert_jitter < 0:
        raise ValueError("--insert-jitter must be zero or greater")
    if args.insert_size - args.insert_jitter < args.read_length:
        raise ValueError(
            "--insert-size minus --insert-jitter must be at least --read-length"
        )
    return primer


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_manifest(
    path: Path,
    *,
    reference: Path,
    records: Sequence[SeqRecord],
    candidates: Sequence[CandidateSite],
    counts: Counter[CandidateSite],
    generated_read_pairs: int,
    args: argparse.Namespace,
    primer: str,
    output_files: Sequence[Path],
) -> None:
    manifest = {
        "format_version": 1,
        "reference": str(reference.resolve()),
        "reference_sha256": sha256(reference),
        "contigs": [
            {"id": record.id, "length_bp": len(record.seq)} for record in records
        ],
        "candidate_model": "Himar1/mariner TA",
        "coordinate_system": "1-based position of the T in each TA dinucleotide",
        "candidate_sites": len(candidates),
        "positive_sites": sum(count > 0 for count in counts.values()),
        "requested_read_pairs": args.read_pairs,
        "generated_read_pairs": generated_read_pairs,
        "read_length": args.read_length,
        "insert_size": args.insert_size,
        "insert_jitter": args.insert_jitter,
        "r1_terminal_primer": primer,
        "seed": args.seed,
        "files": {
            output.name: {"sha256": sha256(output), "bytes": output.stat().st_size}
            for output in output_files
        },
        "warning": (
            "Synthetic smoke-test data only; validate the real transposon primer "
            "and library architecture before biological interpretation."
        ),
    }
    path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    primer = validate_args(args)
    records = load_records(args.reference)
    candidates = enumerate_candidates(records)
    rng = random.Random(args.seed)
    minimum_fragment = args.insert_size - args.insert_jitter
    counts = allocate_counts(
        candidates,
        records,
        read_pairs=args.read_pairs,
        minimum_fragment_length=minimum_fragment,
        rng=rng,
    )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = args.output_dir / "final-dataset.csv"
    read1_path = args.output_dir / "mock_reads1.fastq.gz"
    read2_path = args.output_dir / "mock_reads2.fastq.gz"
    manifest_path = args.output_dir / "mock_manifest.json"

    write_final_dataset(csv_path, candidates, counts)
    generated_read_pairs = write_fastq_pair(
        read1_path,
        read2_path,
        candidates,
        counts,
        records,
        primer=primer,
        read_length=args.read_length,
        insert_size=args.insert_size,
        insert_jitter=args.insert_jitter,
        rng=rng,
    )
    if generated_read_pairs != args.read_pairs:
        raise AssertionError(
            f"Generated {generated_read_pairs} pairs, expected {args.read_pairs}"
        )
    if sum(counts.values()) != generated_read_pairs:
        raise AssertionError("CSV read counts do not equal the FASTQ pair count")

    write_manifest(
        manifest_path,
        reference=args.reference,
        records=records,
        candidates=candidates,
        counts=counts,
        generated_read_pairs=generated_read_pairs,
        args=args,
        primer=primer,
        output_files=(csv_path, read1_path, read2_path),
    )
    print(
        f"Wrote {len(candidates)} TA sites, "
        f"{sum(count > 0 for count in counts.values())} positive sites, and "
        f"{generated_read_pairs} read pairs to {args.output_dir}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
