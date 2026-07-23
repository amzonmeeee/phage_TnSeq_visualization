"""Small, testable building blocks for the optional Tn-Seq processing pipeline.

The visualization package deliberately does not bundle FastQC, fastp, TRANSIT,
BWA, or SeqKit.  They are command-line programs installed by the user (normally
in a conda environment).  This module therefore has two jobs:

* construct reproducible *argument lists* for those programs; and
* turn TRANSIT TPP ``.wig`` output into ordinary insertion-count records.

Command construction itself has no side effects.  The explicit ``run_pipeline``
entry point executes a selected plan, while the lower-level builders remain safe
to test without bioinformatics binaries and make every invocation inspectable in
a run manifest.

TPP already invokes BWA to map the genomic part of a transposon-junction read.
Consequently, ``build_tpp_command`` passes the BWA executable to TPP rather than
constructing a second, duplicate BWA mapping command.
"""

from __future__ import annotations

from dataclasses import dataclass
import gzip
from pathlib import Path
from statistics import pstdev
import subprocess
from typing import Callable, Iterable, Mapping, Sequence


class PipelineError(ValueError):
    """Raised for an invalid processing configuration or malformed count file."""


@dataclass(frozen=True)
class ToolPaths:
    """Names or absolute paths for external executables.

    ``tpp`` may be either an executable (``("tpp.py",)``) or an interpreter plus
    script (``("python", "/path/to/tpp.py")``).  The latter is useful for TRANSIT
    installations that do not put ``tpp.py`` on ``PATH``.
    """

    fastqc: str = "fastqc"
    fastp: str = "fastp"
    tpp: tuple[str, ...] = ("tpp.py",)
    bwa: str = "bwa"
    seqkit: str = "seqkit"


@dataclass(frozen=True)
class ReadPair:
    """One single-end or paired-end FASTQ input, with a stable sample identifier."""

    read1: Path
    read2: Path | None = None
    sample_id: str = "sample"

    def __post_init__(self) -> None:
        if not self.sample_id.strip():
            raise PipelineError("sample_id must not be empty")


@dataclass(frozen=True)
class InsertionCount:
    """Read/template count at one 1-based reference coordinate."""

    contig: str
    position: int
    count: float

    def __post_init__(self) -> None:
        if not self.contig:
            raise PipelineError("contig must not be empty")
        if self.position < 1:
            raise PipelineError("insertion positions must be 1-based positive integers")
        if self.count < 0:
            raise PipelineError("insertion counts must not be negative")

    @property
    def read_count(self) -> float:
        """Semantic alias used by CSV and essentiality consumers."""
        return self.count


@dataclass(frozen=True)
class AveragedInsertionCount:
    """A count averaged over independently seeded, depth-matched subsamples."""

    contig: str
    position: int
    mean_count: float
    sd_count: float
    n_subsamples: int

    @property
    def count(self) -> float:
        """Compatibility alias for consumers that only need the final count."""
        return self.mean_count


@dataclass(frozen=True)
class SubsampleJob:
    """The two SeqKit commands and output pair for one random subsample."""

    replicate: int
    seed: int
    reads: ReadPair
    commands: tuple[tuple[str, ...], ...]


CommandRunner = Callable[[tuple[str, ...]], object]


@dataclass(frozen=True)
class PipelineConfig:
    """Configuration for :func:`run_pipeline`.

    The processing stages are intentionally opt-in/out individually.  Skipping
    FastQC or fastp passes the previous FASTQ files to the following stage.
    Skipping TPP is useful when a caller only wants the generated reference FASTA
    and/or command plan; it cannot produce an insertion count table by itself.
    """

    output_dir: Path
    threads: int = 1
    run_fastqc: bool = True
    run_fastp: bool = True
    run_tpp: bool = True
    minimum_phred: int = 20
    minimum_read_length: int = 20
    adapter_sequence: str | None = None
    adapter_sequence_r2: str | None = None
    tpp_primer: str | None = None
    tpp_mismatches: int = 1
    tpp_mode: str = "himar1"
    tpp_replicon_ids: tuple[str, ...] | None = None
    minimum_mapped_reads: float = 0
    include_unobserved_sites: bool = True
    subsample_depth: int | None = None
    subsample_replicates: int = 1
    subsample_seed: int = 11

    def __post_init__(self) -> None:
        object.__setattr__(self, "output_dir", Path(self.output_dir))
        _require_positive("threads", self.threads)
        _require_phred(self.minimum_phred)
        _require_positive("minimum_read_length", self.minimum_read_length)
        if self.tpp_mismatches < 0:
            raise PipelineError("tpp_mismatches must be zero or greater")
        if self.tpp_mode.lower() not in {"himar1", "tn5"}:
            raise PipelineError("tpp_mode must be 'himar1' or 'tn5'")
        if self.minimum_mapped_reads < 0:
            raise PipelineError("minimum_mapped_reads must be zero or greater")
        if self.subsample_depth is None:
            if self.subsample_replicates != 1:
                raise PipelineError(
                    "subsample_replicates must be 1 when subsample_depth is not set"
                )
        else:
            _require_positive("subsample_depth", self.subsample_depth)
            _require_positive("subsample_replicates", self.subsample_replicates)


@dataclass(frozen=True)
class PipelineResult:
    """Artifacts and parsed counts produced by :func:`run_pipeline`.

    ``counts`` is populated for a non-subsampled TPP run.  When depth matching is
    requested, ``averaged_counts`` contains the final mean/SD records instead.
    Both stay empty when ``PipelineConfig.run_tpp`` is false.
    """

    reference_fasta: Path
    processed_reads: ReadPair
    commands: tuple[tuple[str, ...], ...]
    wig_files: tuple[Path, ...]
    counts: tuple[InsertionCount, ...] = ()
    averaged_counts: tuple[AveragedInsertionCount, ...] = ()


def export_reference_fasta(
    reference_genbank: str | Path,
    output_fasta: str | Path,
) -> Path:
    """Export every GenBank record to the FASTA reference required by TPP/BWA.

    The visualization input remains the original annotated GenBank file.  This
    conversion is only an analysis artifact: it keeps the coordinate system and
    contig identifiers while supplying the FASTA format accepted by TPP and BWA.
    """
    from Bio import SeqIO

    reference_genbank = Path(reference_genbank)
    output_fasta = Path(output_fasta)
    records = list(SeqIO.parse(str(reference_genbank), "genbank"))
    if not records:
        raise PipelineError(f"No GenBank records found in {reference_genbank}")
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    written = SeqIO.write(records, str(output_fasta), "fasta")
    if written != len(records):  # pragma: no cover - defensive Biopython guard
        raise PipelineError(f"Could not write every reference contig to {output_fasta}")
    return output_fasta


def reference_contig_ids(
    reference_genbank: str | Path,
    *,
    replicon_ids: Sequence[str] | None = None,
) -> list[str]:
    """Return the contig name to use for each GenBank record, in file order.

    This is the single source of truth for contig identity across the raw-read
    path.  TPP's WIG header cannot supply it: TPP writes ``chrom=`` as the
    reference FASTA's basename truncated at the first dot, which is the same
    string for every replicon and never the accession.  Contig identity is
    therefore taken from the reference itself, and from the ``-replicon-ids``
    that TPP echoes into its per-replicon output *filenames*.
    """
    from Bio import SeqIO

    records = list(SeqIO.parse(str(reference_genbank), "genbank"))
    if not records:
        raise PipelineError(f"No GenBank records found in {reference_genbank}")
    if replicon_ids is None:
        contigs = [record.id for record in records]
    else:
        if len(replicon_ids) != len(records):
            raise PipelineError(
                "replicon_ids must contain exactly one ID for each GenBank record"
            )
        contigs = [str(item).strip() for item in replicon_ids]
    for contig in contigs:
        if not contig:
            raise PipelineError("reference contig has no identifier")
    if len(set(contigs)) != len(contigs):
        raise PipelineError("reference contig identifiers must be unique")
    return contigs


def candidate_insertion_sites(
    reference_genbank: str | Path,
    *,
    mode: str = "himar1",
    replicon_ids: Sequence[str] | None = None,
) -> dict[str, list[int]]:
    """Return every candidate coordinate at which a library could be observed.

    Sassetti/Himar1 libraries are constrained to TA dinucleotides, whereas Tn5
    is treated as unconstrained and therefore includes every reference base.  The
    result is designed for :func:`fill_missing_candidate_sites`, ensuring that a
    zero-count candidate stays explicit for downstream essentiality classification.

    When TPP is given ``-replicon-ids``, use the same ordered IDs here so WIG
    contigs and reference candidates share the same names.
    """
    from Bio import SeqIO

    mode = mode.lower()
    if mode not in {"himar1", "tn5"}:
        raise PipelineError("mode must be 'himar1' or 'tn5'")
    records = list(SeqIO.parse(str(reference_genbank), "genbank"))
    if not records:
        raise PipelineError(f"No GenBank records found in {reference_genbank}")
    contigs = reference_contig_ids(reference_genbank, replicon_ids=replicon_ids)

    candidates: dict[str, list[int]] = {}
    for contig, record in zip(contigs, records):
        sequence = str(record.seq)
        if mode == "himar1":
            # TA is palindromic, so one scan fully enumerates insertion positions.
            from .transposon import find_insertion_sites

            candidates[contig] = find_insertion_sites(sequence, "TA", both_strands=False)
        else:
            candidates[contig] = list(range(1, len(sequence) + 1))
    return candidates


def run_pipeline(
    reference_genbank: str | Path,
    reads: ReadPair,
    config: PipelineConfig,
    *,
    tools: ToolPaths = ToolPaths(),
    runner: CommandRunner | None = None,
) -> PipelineResult:
    """Run the selected raw-read processing stages and parse the resulting WIG file(s).

    ``runner`` receives an immutable argv tuple for each command.  It exists both
    for application-specific progress/log handling and for dry unit tests; when
    omitted, commands are executed with :func:`subprocess.run` and ``check=True``.
    Commands are never sent through a shell.

    TPP is the only mapping stage: it receives ``-bwa`` and carries out the BWA
    mapping internally.  If a depth match is requested, every SeqKit subsample is
    processed through TPP independently, thresholded, and then averaged.
    """
    output_dir = config.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    reference_fasta = export_reference_fasta(reference_genbank, output_dir / "reference.fasta")
    execute = runner or _default_runner
    commands: list[tuple[str, ...]] = []

    def run(command: tuple[str, ...]) -> None:
        commands.append(command)
        execute(command)

    current_reads = reads
    if config.run_fastqc:
        fastqc_dir = output_dir / "01_fastqc"
        fastqc_dir.mkdir(parents=True, exist_ok=True)
        run(build_fastqc_command(current_reads, fastqc_dir, tools=tools, threads=config.threads))

    if config.run_fastp:
        fastp_dir = output_dir / "02_fastp"
        fastp_dir.mkdir(parents=True, exist_ok=True)
        fastp_input = current_reads
        current_reads = _trimmed_read_pair(fastp_input, fastp_dir)
        run(
            build_fastp_command(
                fastp_input,
                current_reads,
                html_report=fastp_dir / f"{reads.sample_id}.fastp.html",
                json_report=fastp_dir / f"{reads.sample_id}.fastp.json",
                minimum_phred=config.minimum_phred,
                minimum_length=config.minimum_read_length,
                threads=config.threads,
                adapter_sequence=config.adapter_sequence,
                adapter_sequence_r2=config.adapter_sequence_r2,
                tools=tools,
            )
        )

    if not config.run_tpp:
        return PipelineResult(
            reference_fasta=reference_fasta,
            processed_reads=current_reads,
            commands=tuple(commands),
            wig_files=(),
        )

    contigs = reference_contig_ids(reference_genbank, replicon_ids=config.tpp_replicon_ids)
    if len(contigs) > 1 and config.tpp_replicon_ids is None:
        raise PipelineError(
            f"{reference_genbank} contains {len(contigs)} records, so TPP needs "
            "--tpp-replicon-ids with one ID per record; without it TPP cannot write "
            "separate per-contig WIG files and the counts cannot be told apart"
        )

    candidate_sites = (
        candidate_insertion_sites(
            reference_genbank,
            mode=config.tpp_mode,
            replicon_ids=config.tpp_replicon_ids,
        )
        if config.include_unobserved_sites
        else None
    )
    tpp_dir = output_dir / "03_tpp"
    tpp_dir.mkdir(parents=True, exist_ok=True)
    if config.subsample_depth is None:
        prefix = tpp_dir / reads.sample_id
        run(
            build_tpp_command(
                reference_fasta,
                current_reads,
                prefix,
                primer=config.tpp_primer,
                mismatches=config.tpp_mismatches,
                mode=config.tpp_mode,
                replicon_ids=config.tpp_replicon_ids,
                tools=tools,
            )
        )
        wig_pairs = _find_tpp_wig_files(prefix, contigs)
        rows = _parse_wig_files(wig_pairs)
        if candidate_sites is not None:
            rows = fill_missing_candidate_sites(rows, candidate_sites)
        rows = apply_read_count_threshold(rows, config.minimum_mapped_reads)
        return PipelineResult(
            reference_fasta=reference_fasta,
            processed_reads=current_reads,
            commands=tuple(commands),
            wig_files=tuple(path for path, _contig in wig_pairs),
            counts=tuple(rows),
        )

    subsample_dir = output_dir / "03_subsamples"
    subsample_dir.mkdir(parents=True, exist_ok=True)
    jobs = build_subsample_jobs(
        current_reads,
        subsample_dir,
        depth=config.subsample_depth,
        replicates=config.subsample_replicates,
        seed=config.subsample_seed,
        tools=tools,
    )
    all_wig_files: list[Path] = []
    replicate_tables: list[list[InsertionCount]] = []
    for job in jobs:
        for command in job.commands:
            run(command)
        prefix = tpp_dir / f"{reads.sample_id}.subsample-{job.replicate:03d}"
        run(
            build_tpp_command(
                reference_fasta,
                job.reads,
                prefix,
                primer=config.tpp_primer,
                mismatches=config.tpp_mismatches,
                mode=config.tpp_mode,
                replicon_ids=config.tpp_replicon_ids,
                tools=tools,
            )
        )
        wig_pairs = _find_tpp_wig_files(prefix, contigs)
        all_wig_files.extend(path for path, _contig in wig_pairs)
        rows = _parse_wig_files(wig_pairs)
        if candidate_sites is not None:
            rows = fill_missing_candidate_sites(rows, candidate_sites)
        replicate_tables.append(rows)

    averaged = average_subsampled_counts(
        replicate_tables,
        minimum_count=config.minimum_mapped_reads,
    )
    return PipelineResult(
        reference_fasta=reference_fasta,
        processed_reads=current_reads,
        commands=tuple(commands),
        wig_files=tuple(all_wig_files),
        averaged_counts=tuple(averaged),
    )


def build_fastqc_command(
    reads: ReadPair,
    output_dir: str | Path,
    *,
    tools: ToolPaths = ToolPaths(),
    threads: int = 1,
) -> tuple[str, ...]:
    """Return the FastQC command for the input FASTQ file(s)."""
    _require_positive("threads", threads)
    cmd = [tools.fastqc, "--outdir", str(Path(output_dir)), "--threads", str(threads)]
    cmd.append(str(reads.read1))
    if reads.read2 is not None:
        cmd.append(str(reads.read2))
    return tuple(cmd)


def build_fastp_command(
    reads: ReadPair,
    output_reads: ReadPair,
    *,
    html_report: str | Path,
    json_report: str | Path,
    minimum_phred: int = 20,
    minimum_length: int = 20,
    threads: int = 1,
    adapter_sequence: str | None = None,
    adapter_sequence_r2: str | None = None,
    tools: ToolPaths = ToolPaths(),
) -> tuple[str, ...]:
    """Return a fastp adapter-trimming and quality-filtering command.

    fastp enables adapter trimming by default.  For paired-end data we also add
    ``--detect_adapter_for_pe`` so it attempts overlap-based adapter detection;
    callers can supply known adapter sequences to override that discovery.
    """
    _require_phred(minimum_phred)
    _require_positive("minimum_length", minimum_length)
    _require_positive("threads", threads)
    _validate_pair_shape(reads, output_reads)
    if adapter_sequence_r2 and reads.read2 is None:
        raise PipelineError("adapter_sequence_r2 requires paired-end reads")

    cmd = [
        tools.fastp,
        "--in1", str(reads.read1),
        "--out1", str(output_reads.read1),
        "--qualified_quality_phred", str(minimum_phred),
        "--length_required", str(minimum_length),
        "--thread", str(threads),
        "--html", str(Path(html_report)),
        "--json", str(Path(json_report)),
    ]
    if reads.read2 is not None:
        # _validate_pair_shape guarantees output_reads.read2 is also present.
        cmd.extend([
            "--in2", str(reads.read2),
            "--out2", str(output_reads.read2),
            "--detect_adapter_for_pe",
        ])
    if adapter_sequence:
        cmd.extend(["--adapter_sequence", adapter_sequence])
    if adapter_sequence_r2:
        cmd.extend(["--adapter_sequence_r2", adapter_sequence_r2])
    return tuple(cmd)


def build_tpp_command(
    reference_fasta: str | Path,
    reads: ReadPair,
    output_prefix: str | Path,
    *,
    primer: str | None = None,
    mismatches: int = 1,
    mode: str = "himar1",
    replicon_ids: Sequence[str] | None = None,
    tools: ToolPaths = ToolPaths(),
) -> tuple[str, ...]:
    """Return a TRANSIT TPP command that performs its integrated BWA mapping.

    ``mode`` is intentionally limited to the two library families exposed by
    this application: ``"himar1"`` (TRANSIT v3 protocol ``Sassetti``, TA-site
    output) and ``"tn5"`` (protocol ``Tn5``, all genomic positions).  A
    non-standard IR--cargo--IR construct must still provide the actual terminal
    sequence using ``primer``; this helper does not pretend to infer it from the
    cargo sequence.
    """
    if mismatches < 0:
        raise PipelineError("mismatches must be zero or greater")
    mode = mode.lower()
    if mode not in {"himar1", "tn5"}:
        raise PipelineError("mode must be 'himar1' or 'tn5'")
    if primer is not None and not primer.strip():
        raise PipelineError("primer must not be blank when supplied")
    if not tools.tpp:
        raise PipelineError("at least one TPP command token is required")

    protocol = {"himar1": "Sassetti", "tn5": "Tn5"}[mode]
    cmd = [
        *tools.tpp,
        "-bwa", tools.bwa,
        "-ref", str(Path(reference_fasta)),
        "-reads1", str(reads.read1),
    ]
    if reads.read2 is not None:
        cmd.extend(["-reads2", str(reads.read2)])
    cmd.extend([
        "-output", str(Path(output_prefix)),
        "-mismatches", str(mismatches),
        "-protocol", protocol,
    ])
    if primer is not None:
        cmd.extend(["-primer", primer])
    if replicon_ids:
        if any(not item.strip() for item in replicon_ids):
            raise PipelineError("replicon_ids must not contain blank IDs")
        cmd.extend(["-replicon-ids", ",".join(replicon_ids)])
    return tuple(cmd)


def build_subsample_jobs(
    reads: ReadPair,
    output_dir: str | Path,
    *,
    depth: int,
    replicates: int,
    seed: int = 11,
    tools: ToolPaths = ToolPaths(),
) -> list[SubsampleJob]:
    """Plan reproducible, depth-matched SeqKit ``sample2`` jobs.

    Paired reads are sampled in separate commands using the *same* seed.  SeqKit
    documents that this produces the same selected read pairs in each FASTQ;
    `sample2 -2` uses two passes to keep memory bounded on deep sequencing data.
    No command is run by this function.
    """
    _require_positive("depth", depth)
    _require_positive("replicates", replicates)
    output_dir = Path(output_dir)
    jobs: list[SubsampleJob] = []
    for replicate in range(1, replicates + 1):
        replicate_seed = seed + replicate - 1
        out1 = output_dir / _subsample_name(reads.sample_id, replicate, "R1", reads.read1)
        out2 = (
            output_dir / _subsample_name(reads.sample_id, replicate, "R2", reads.read2)
            if reads.read2 is not None
            else None
        )
        out_reads = ReadPair(out1, out2, reads.sample_id)
        commands = [
            _build_seqkit_sample_command(
                reads.read1, out1, depth, replicate_seed, tools.seqkit
            )
        ]
        if reads.read2 is not None:
            commands.append(
                _build_seqkit_sample_command(
                    reads.read2, out2, depth, replicate_seed, tools.seqkit
                )
            )
        jobs.append(
            SubsampleJob(
                replicate=replicate,
                seed=replicate_seed,
                reads=out_reads,
                commands=tuple(commands),
            )
        )
    return jobs


def parse_tpp_wig(
    path: str | Path,
    *,
    default_contig: str | None = None,
    contig: str | None = None,
) -> list[InsertionCount]:
    """Parse a TPP WIG file into sorted, coalesced insertion-count records.

    Both ``variableStep`` (the usual TPP form) and ``fixedStep`` WIG sections are
    accepted.  Duplicate coordinate lines are summed defensively.  A headerless
    two-column file is also accepted when a single ``default_contig`` is supplied.

    ``contig`` overrides whatever the header declares.  TPP output requires it:
    TPP writes the reference FASTA's basename into ``chrom=`` rather than a
    contig name, so every replicon of a multi-contig run claims the same one.
    """
    path = Path(path)
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        return parse_tpp_wig_lines(handle, default_contig=default_contig, contig=contig)


def parse_tpp_wig_lines(
    lines: Iterable[str],
    *,
    default_contig: str | None = None,
    contig: str | None = None,
) -> list[InsertionCount]:
    """Line-oriented implementation of :func:`parse_tpp_wig` for unit tests."""
    if contig is not None and not contig.strip():
        raise PipelineError("contig override must not be blank")
    override = contig.strip() if contig is not None else None
    contig = override or default_contig
    fixed_position: int | None = None
    fixed_step = 1
    counts: dict[tuple[str, int], float] = {}

    for line_number, raw_line in enumerate(lines, start=1):
        line = raw_line.strip()
        if not line or line.startswith("#") or line.startswith("track"):
            continue
        fields = line.split()
        directive = fields[0].lower()
        if directive in {"variablestep", "fixedstep"}:
            attrs = _wig_attributes(fields[1:], line_number)
            declared = attrs.get("chrom") or attrs.get("chromosome")
            if not declared:
                raise PipelineError(f"WIG line {line_number}: step declaration lacks chrom=")
            contig = override or declared
            if directive == "fixedstep":
                try:
                    fixed_position = int(attrs.get("start", "1"))
                    fixed_step = int(attrs.get("step", "1"))
                except ValueError as exc:
                    raise PipelineError(
                        f"WIG line {line_number}: fixedStep start and step must be integers"
                    ) from exc
                if fixed_position < 1 or fixed_step < 1:
                    raise PipelineError(
                        f"WIG line {line_number}: fixedStep start and step must be positive"
                    )
            else:
                fixed_position = None
            continue

        if len(fields) == 2:
            try:
                position = int(fields[0])
                value = float(fields[1])
            except ValueError as exc:
                raise PipelineError(
                    f"WIG line {line_number}: expected '<position> <count>'"
                ) from exc
        elif len(fields) == 1 and fixed_position is not None:
            position = fixed_position
            try:
                value = float(fields[0])
            except ValueError as exc:
                raise PipelineError(f"WIG line {line_number}: invalid count") from exc
            fixed_position += fixed_step
        else:
            raise PipelineError(
                f"WIG line {line_number}: expected a WIG step declaration or count value"
            )

        if not contig:
            raise PipelineError(
                f"WIG line {line_number}: no contig is known; add a step header or default_contig"
            )
        if position < 1:
            raise PipelineError(f"WIG line {line_number}: position must be positive")
        if value < 0:
            raise PipelineError(f"WIG line {line_number}: count must not be negative")
        counts[(contig, position)] = counts.get((contig, position), 0.0) + value

    return [
        InsertionCount(contig, position, count)
        for (contig, position), count in sorted(counts.items())
    ]


def fill_missing_candidate_sites(
    counts: Iterable[InsertionCount],
    candidate_sites: Mapping[str, Iterable[int]],
) -> list[InsertionCount]:
    """Return every expected candidate site, filling unobserved sites with zero.

    This is most useful for a Himar1/TA library: scan the reference for TA sites,
    then left-join the TPP counts so genuinely unhit sites remain visible rather
    than disappearing from the final dataset.
    """
    count_map = _count_map(counts)
    filled: list[InsertionCount] = []
    for contig, positions in candidate_sites.items():
        if not contig:
            raise PipelineError("candidate-site contig must not be empty")
        for position in sorted(set(positions)):
            if position < 1:
                raise PipelineError("candidate-site positions must be positive")
            filled.append(InsertionCount(contig, position, count_map.get((contig, position), 0.0)))
    return sorted(filled, key=lambda row: (row.contig, row.position))


def apply_read_count_threshold(
    counts: Iterable[InsertionCount],
    minimum_count: float,
) -> list[InsertionCount]:
    """Set counts strictly below ``minimum_count`` to zero.

    The strict comparison implements the requested "less than threshold" rule:
    a site whose count equals the threshold remains non-zero.
    """
    if minimum_count < 0:
        raise PipelineError("minimum_count must be zero or greater")
    return [
        InsertionCount(row.contig, row.position, 0.0 if row.count < minimum_count else row.count)
        for row in counts
    ]


def average_subsampled_counts(
    replicate_counts: Sequence[Iterable[InsertionCount]],
    *,
    minimum_count: float = 0,
) -> list[AveragedInsertionCount]:
    """Threshold each depth-matched replicate, then calculate mean and population SD.

    Missing coordinates in a replicate are interpreted as zero counts.  Applying
    the threshold *before* the average preserves the meaning of the chosen cutoff
    for each independently resampled library.
    """
    if not replicate_counts:
        raise PipelineError("at least one subsample count table is required")
    if minimum_count < 0:
        raise PipelineError("minimum_count must be zero or greater")

    maps = [_count_map(apply_read_count_threshold(rows, minimum_count)) for rows in replicate_counts]
    sites = sorted({site for table in maps for site in table})
    n = len(maps)
    averaged: list[AveragedInsertionCount] = []
    for contig, position in sites:
        values = [table.get((contig, position), 0.0) for table in maps]
        averaged.append(
            AveragedInsertionCount(
                contig=contig,
                position=position,
                mean_count=sum(values) / n,
                sd_count=pstdev(values) if n > 1 else 0.0,
                n_subsamples=n,
            )
        )
    return averaged


def parse_wig(path: str | Path, *, default_contig: str | None = None) -> list[InsertionCount]:
    """Plain-name alias for :func:`parse_tpp_wig` used by CLI integrations."""
    return parse_tpp_wig(path, default_contig=default_contig)


def threshold_counts(
    counts: Iterable[InsertionCount],
    minimum_count: float,
) -> list[InsertionCount]:
    """Plain-name alias for :func:`apply_read_count_threshold`."""
    return apply_read_count_threshold(counts, minimum_count)


def _require_positive(name: str, value: int) -> None:
    if value < 1:
        raise PipelineError(f"{name} must be greater than zero")


def _require_phred(value: int) -> None:
    if not 0 <= value <= 93:
        raise PipelineError("minimum_phred must be between 0 and 93")


def _validate_pair_shape(reads: ReadPair, output_reads: ReadPair) -> None:
    if (reads.read2 is None) != (output_reads.read2 is None):
        raise PipelineError("input and output reads must both be single-end or both be paired-end")


def _wig_attributes(fields: Sequence[str], line_number: int) -> dict[str, str]:
    """Parse ``key=value`` step attributes, tolerating TPP's comma separators.

    A multi-replicon TPP run writes ``variableStep chrom=<ref>, replicon=<i>``.
    The trailing comma is not valid WIG, so it is stripped rather than becoming
    part of the contig name.
    """
    attrs: dict[str, str] = {}
    for field in fields:
        field = field.rstrip(",")
        if not field:
            continue
        if "=" not in field:
            raise PipelineError(f"WIG line {line_number}: malformed step attribute {field!r}")
        key, value = field.split("=", 1)
        attrs[key.lower()] = value.strip()
    return attrs


def _count_map(counts: Iterable[InsertionCount]) -> dict[tuple[str, int], float]:
    result: dict[tuple[str, int], float] = {}
    for row in counts:
        key = (row.contig, row.position)
        result[key] = result.get(key, 0.0) + row.count
    return result


def _trimmed_read_pair(reads: ReadPair, output_dir: Path) -> ReadPair:
    """Choose deterministic fastp output names while preserving FASTQ compression."""
    out1 = output_dir / f"{reads.sample_id}.trimmed.R1{_fastq_suffix(reads.read1)}"
    out2 = (
        output_dir / f"{reads.sample_id}.trimmed.R2{_fastq_suffix(reads.read2)}"
        if reads.read2 is not None
        else None
    )
    return ReadPair(out1, out2, reads.sample_id)


def _find_tpp_wig_files(prefix: Path, contigs: Sequence[str]) -> list[tuple[Path, str]]:
    """Pair each TPP WIG with the contig it holds, using TPP's own naming rule.

    TPP writes ``<prefix>.wig`` for a single-replicon reference and
    ``<prefix>_<replicon_id>.wig`` for each replicon otherwise, where the IDs are
    exactly what was passed to ``-replicon-ids``.  The filename is the only place
    TPP records which replicon a file belongs to; its ``chrom=`` header does not.
    """
    if not contigs:
        raise PipelineError("at least one reference contig is required")

    if len(contigs) == 1:
        path = prefix.with_name(f"{prefix.name}.wig")
        if not path.exists():
            raise PipelineError(
                f"TPP completed but {path.name!r} was not found in {prefix.parent}"
            )
        return [(path, contigs[0])]

    pairs: list[tuple[Path, str]] = []
    missing: list[str] = []
    for contig in contigs:
        path = prefix.with_name(f"{prefix.name}_{contig}.wig")
        if path.exists():
            pairs.append((path, contig))
        else:
            missing.append(path.name)
    if missing:
        raise PipelineError(
            "TPP completed but these per-replicon WIG files were not found in "
            f"{prefix.parent}: {', '.join(missing)}. Check that --tpp-replicon-ids "
            "matches the reference record order."
        )
    return pairs


def _parse_wig_files(pairs: Iterable[tuple[Path, str]]) -> list[InsertionCount]:
    """Parse per-replicon WIGs, forcing each one's contig name."""
    rows: list[InsertionCount] = []
    for path, contig in pairs:
        rows.extend(parse_tpp_wig(path, contig=contig))
    count_map = _count_map(rows)
    return [
        InsertionCount(contig, position, count)
        for (contig, position), count in sorted(count_map.items())
    ]


def _default_runner(command: tuple[str, ...]) -> None:
    """Execute one external command without a shell so paths remain safely quoted."""
    try:
        subprocess.run(command, check=True)
    except FileNotFoundError as exc:
        raise PipelineError(
            f"Required external program was not found: {command[0]!r}. "
            "Install it or configure its path explicitly."
        ) from exc
    except subprocess.CalledProcessError as exc:
        raise PipelineError(
            f"External pipeline command failed with exit code {exc.returncode}: {' '.join(command)}"
        ) from exc


def _subsample_name(sample_id: str, replicate: int, read_label: str, source: Path | None) -> str:
    suffix = _fastq_suffix(source) if source is not None else ".fastq.gz"
    return f"{sample_id}.subsample-{replicate:03d}.{read_label}{suffix}"


def _fastq_suffix(path: Path) -> str:
    name = path.name.lower()
    for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if name.endswith(suffix):
            return suffix
    # Preserve an unfamiliar extension rather than silently changing compression.
    return path.suffix or ".fastq.gz"


def _build_seqkit_sample_command(
    source: Path,
    output: Path | None,
    depth: int,
    seed: int,
    seqkit: str,
) -> tuple[str, ...]:
    if output is None:  # Defensive; only paired jobs call this for R2.
        raise PipelineError("subsample output path must not be None")
    return (
        seqkit,
        "sample2",
        "-n", str(depth),
        "-2",
        "-s", str(seed),
        "-o", str(output),
        str(source),
    )
