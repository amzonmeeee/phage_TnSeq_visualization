"""Tests for optional external-tool pipeline planning and TPP count parsing.

These tests only inspect generated argv lists and small in-memory WIG snippets.
They deliberately do not require FastQC, fastp, TRANSIT, BWA, or SeqKit.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from phage_tnseq_viz.pipeline import (
    InsertionCount,
    PipelineConfig,
    PipelineError,
    ReadPair,
    ToolPaths,
    apply_read_count_threshold,
    average_subsampled_counts,
    build_fastp_command,
    build_subsample_jobs,
    build_tpp_command,
    candidate_insertion_sites,
    fill_missing_candidate_sites,
    parse_tpp_wig_lines,
    run_pipeline,
)


EXAMPLES = Path(__file__).resolve().parent.parent / "examples"


def test_fastp_command_contains_requested_quality_and_adapter_options():
    reads = ReadPair(Path("raw_R1.fastq.gz"), Path("raw_R2.fastq.gz"), "phage-a")
    trimmed = ReadPair(Path("trim_R1.fastq.gz"), Path("trim_R2.fastq.gz"), "phage-a")

    cmd = build_fastp_command(
        reads,
        trimmed,
        html_report=Path("fastp.html"),
        json_report=Path("fastp.json"),
        minimum_phred=28,
        minimum_length=45,
        threads=6,
        adapter_sequence="AGATCGGAAGAGC",
        adapter_sequence_r2="AGATCGGAAGAGC",
        tools=ToolPaths(fastp="/tools/fastp"),
    )

    assert cmd[:5] == ("/tools/fastp", "--in1", "raw_R1.fastq.gz", "--out1", "trim_R1.fastq.gz")
    assert ("--qualified_quality_phred", "28") == _pair_after(cmd, "--qualified_quality_phred")
    assert ("--length_required", "45") == _pair_after(cmd, "--length_required")
    assert ("--thread", "6") == _pair_after(cmd, "--thread")
    assert "--detect_adapter_for_pe" in cmd
    assert ("--adapter_sequence", "AGATCGGAAGAGC") == _pair_after(cmd, "--adapter_sequence")
    assert ("--adapter_sequence_r2", "AGATCGGAAGAGC") == _pair_after(cmd, "--adapter_sequence_r2")


def test_tpp_command_delegates_bwa_to_tpp_and_keeps_custom_prefix():
    reads = ReadPair(Path("trim_R1.fastq.gz"), Path("trim_R2.fastq.gz"), "phage-a")
    tools = ToolPaths(tpp=("python", "/opt/transit/tpp.py"), bwa="/opt/bwa")

    cmd = build_tpp_command(
        Path("reference.fasta"),
        reads,
        Path("out/phage-a"),
        primer="ACGTACGT",
        mismatches=2,
        mode="tn5",
        replicon_ids=("contig_1", "contig_2"),
        tools=tools,
    )

    assert cmd[:2] == ("python", "/opt/transit/tpp.py")
    assert ("-bwa", "/opt/bwa") == _pair_after(cmd, "-bwa")
    assert ("-ref", "reference.fasta") == _pair_after(cmd, "-ref")
    assert ("-reads1", "trim_R1.fastq.gz") == _pair_after(cmd, "-reads1")
    assert ("-reads2", "trim_R2.fastq.gz") == _pair_after(cmd, "-reads2")
    assert ("-primer", "ACGTACGT") == _pair_after(cmd, "-primer")
    assert ("-mismatches", "2") == _pair_after(cmd, "-mismatches")
    assert ("-replicon-ids", "contig_1,contig_2") == _pair_after(cmd, "-replicon-ids")
    assert ("-protocol", "Tn5") == _pair_after(cmd, "-protocol")
    assert "mem" not in cmd  # no duplicate standalone BWA-MEM stage


def test_tpp_rejects_unsupported_library_mode():
    with pytest.raises(PipelineError, match="himar1"):
        build_tpp_command("ref.fa", ReadPair(Path("r1.fq")), "out", mode="custom")


def test_seqkit_subsampling_uses_matching_seed_for_each_paired_fastq():
    reads = ReadPair(Path("trim_R1.fastq.gz"), Path("trim_R2.fastq.gz"), "phage-a")

    jobs = build_subsample_jobs(reads, Path("subsamples"), depth=5000, replicates=2, seed=101)

    assert [job.seed for job in jobs] == [101, 102]
    assert [job.replicate for job in jobs] == [1, 2]
    assert jobs[0].reads == ReadPair(
        Path("subsamples/phage-a.subsample-001.R1.fastq.gz"),
        Path("subsamples/phage-a.subsample-001.R2.fastq.gz"),
        "phage-a",
    )
    first_r1, first_r2 = jobs[0].commands
    assert first_r1[:2] == ("seqkit", "sample2")
    assert first_r2[:2] == ("seqkit", "sample2")
    assert ("-n", "5000") == _pair_after(first_r1, "-n")
    assert "-2" in first_r1  # SeqKit two-pass mode for deep FASTQ
    assert ("-s", "101") == _pair_after(first_r1, "-s")
    assert ("-s", "101") == _pair_after(first_r2, "-s")
    assert first_r1[-1] == "trim_R1.fastq.gz"
    assert first_r2[-1] == "trim_R2.fastq.gz"


def test_parse_tpp_variable_step_wig_sorts_and_coalesces_duplicate_sites():
    rows = parse_tpp_wig_lines(
        [
            "track type=wiggle_0\n",
            "variableStep chrom=phage_b\n",
            "8 2\n",
            "4 1\n",
            "4 3\n",
            "variableStep chrom=phage_a\n",
            "2 7\n",
        ]
    )

    assert rows == [
        InsertionCount("phage_a", 2, 7.0),
        InsertionCount("phage_b", 4, 4.0),
        InsertionCount("phage_b", 8, 2.0),
    ]


def test_parse_tpp_fixed_step_wig_and_headerless_default_contig():
    fixed = parse_tpp_wig_lines(["fixedStep chrom=phage start=5 step=4\n", "1\n", "3\n"])
    plain = parse_tpp_wig_lines(["3 9\n"], default_contig="only_contig")

    assert fixed == [InsertionCount("phage", 5, 1.0), InsertionCount("phage", 9, 3.0)]
    assert plain == [InsertionCount("only_contig", 3, 9.0)]


def test_fill_threshold_and_subsample_average_preserve_zeroes():
    candidate_rows = fill_missing_candidate_sites(
        [InsertionCount("phage", 10, 1), InsertionCount("phage", 30, 7)],
        {"phage": [10, 20, 30]},
    )
    thresholded = apply_read_count_threshold(candidate_rows, minimum_count=2)

    assert thresholded == [
        InsertionCount("phage", 10, 0.0),
        InsertionCount("phage", 20, 0.0),
        InsertionCount("phage", 30, 7.0),
    ]

    averaged = average_subsampled_counts(
        [
            [InsertionCount("phage", 10, 1), InsertionCount("phage", 20, 4)],
            [InsertionCount("phage", 10, 4), InsertionCount("phage", 30, 3)],
        ],
        minimum_count=2,
    )

    # The per-replicate threshold turns the first phage:10 count into zero;
    # a site missing from a replicate is also represented as zero.
    assert [(row.position, row.mean_count, row.sd_count, row.n_subsamples) for row in averaged] == [
        (10, 2.0, 2.0, 2),
        (20, 2.0, 2.0, 2),
        (30, 1.5, 1.5, 2),
    ]


def test_run_pipeline_uses_injected_runner_and_parses_its_tpp_output(tmp_path: Path):
    invoked: list[tuple[str, ...]] = []

    def fake_runner(command: tuple[str, ...]) -> None:
        invoked.append(command)
        if "-output" in command:
            prefix = Path(command[command.index("-output") + 1])
            prefix.with_suffix(".wig").write_text(
                "variableStep chrom=NC_TEST01\n10 1\n20 3\n", encoding="utf-8"
            )

    result = run_pipeline(
        EXAMPLES / "example_annotated.gbk",
        ReadPair(Path("raw_R1.fastq.gz"), sample_id="synthetic"),
        PipelineConfig(
            output_dir=tmp_path / "run",
            run_fastqc=False,
            run_fastp=False,
            minimum_mapped_reads=2,
            include_unobserved_sites=False,
        ),
        runner=fake_runner,
    )

    assert result.reference_fasta.exists()
    assert len(invoked) == 1
    assert invoked[0] == result.commands[0]
    assert result.counts == (
        InsertionCount("NC_TEST01", 10, 0.0),
        InsertionCount("NC_TEST01", 20, 3.0),
    )


def test_pipeline_defaults_to_explicit_zeroes_for_unobserved_himar1_ta_sites(tmp_path: Path):
    candidates = candidate_insertion_sites(EXAMPLES / "example_annotated.gbk")
    hit = candidates["NC_TEST01"][0]

    def fake_runner(command: tuple[str, ...]) -> None:
        if "-output" in command:
            prefix = Path(command[command.index("-output") + 1])
            prefix.with_suffix(".wig").write_text(
                f"variableStep chrom=NC_TEST01\n{hit} 3\n", encoding="utf-8"
            )

    result = run_pipeline(
        EXAMPLES / "example_annotated.gbk",
        ReadPair(Path("raw_R1.fastq.gz"), sample_id="synthetic"),
        PipelineConfig(output_dir=tmp_path / "run", run_fastqc=False, run_fastp=False),
        runner=fake_runner,
    )

    assert len(result.counts) == len(candidates["NC_TEST01"])
    assert next(row for row in result.counts if row.position == hit).read_count == 3
    assert any(row.read_count == 0 for row in result.counts)


def _pair_after(command: tuple[str, ...], option: str) -> tuple[str, str]:
    index = command.index(option)
    return command[index], command[index + 1]
