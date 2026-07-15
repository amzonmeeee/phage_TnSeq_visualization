# Phage Tn-Seq Visualization Tool

`phage-tnseq-viz` turns an annotated phage GenBank reference and Tn-Seq data
into a linear genome map, a per-insertion-site CSV, and gene-level
essentiality calls. It can either orchestrate optional raw-read processing or
start directly from a final count table produced elsewhere.

The tool is designed for transposon-enriched Illumina libraries: sheared
genomes, adapter ligation, transposon-specific PCR, and junction sequencing.
It does not attempt to infer a non-standard IR–cargo–IR primer sequence; that
sequence is an experimental input and must be supplied to TPP when needed.

## What is new in v0.2

- Optional FastQC → fastp → TRANSIT TPP + BWA → SeqKit preprocessing path.
- A direct `--final-dataset` path that requires no external bioinformatics
  programs.
- Blue per-site read-count bars above gene arrows.
- Essentiality-coloured gene fills, while PHROG functional categories remain
  arrow borders.
- Normalised site-level and gene-level CSV output, processing logs, and a
  reproducibility manifest.
- A dependency-free implementation of the supplied Harms-lab R classification
  rules, plus a documented custom classifier interface.

## The two input paths

```text
raw R1/R2 FASTQ
  └─ FastQC (optional) → fastp (optional) → TPP + BWA → cutoff
       └─ SeqKit depth-matched subsamples (optional) → final site CSV → map + gene calls

already processed final site CSV
  └──────────────────────────────────────────────────────────────────→ map + gene calls
```

TPP already calls BWA internally to map the genomic portion of junction reads,
tabulate insertion/template counts, and write WIG output. Therefore this tool
uses one integrated **TPP + BWA** stage rather than incorrectly mapping the
same reads a second time with a separate BWA command.

## Installation

```bash
git clone <this-repo>
cd phage_TnSeq_visualization
pip install .
```

For development:

```bash
pip install -e ".[dev]"
pytest
```

The Python package includes Biopython, pyGenomeViz, matplotlib, and
pyrodigal-gv. The raw-read route deliberately does **not** install external
programs for you. Install the tools you plan to run (normally in a dedicated
conda environment) and make them available on `PATH`, or pass their paths with
`--fastqc-bin`, `--fastp-bin`, `--tpp-bin`, `--bwa-bin`, and `--seqkit-bin`.

| Optional stage | External program | What it does here |
|---|---|---|
| Initial QC | FastQC | Generates initial FASTQ quality reports. |
| Trim/filter | fastp | Adapter trimming, qualified-base cutoff, minimum length. |
| Junction calling/counting | TRANSIT TPP + BWA | Finds transposon prefixes, maps genomic junctions, and writes per-site WIG counts. |
| Depth matching | SeqKit `sample2` | Seeded random read subsampling before independent TPP runs. |

## Quick, reproducible demo

The example generator makes an annotated synthetic phage and a **complete**
Himar1/TA final dataset (including zero-count TA sites). It is the fastest way
to test the full visualisation/CSV path without external software.

```bash
python examples/make_example.py --output-dir demo

phage-tnseq-viz plot demo/example_annotated.gbk \
  --final-dataset demo/example_final_sites.csv \
  --output demo/example_map.png \
  --csv-dir demo/results
```

This writes `example_map.png`, `example_map_sites.csv`, and
`example_map_gene_essentiality.csv`.

## Use a final dataset you already processed

If you processed the reads with another validated workflow, do not rerun the
built-in pipeline. Supply the final count table directly:

```bash
phage-tnseq-viz plot phold_annotated_phage.gbk \
  --final-dataset my_final_sites.csv \
  --output phage_tnseq.png \
  --csv-dir results
```

With final data, the plot defaults to the requested blue measured-count bars
and essentiality arrows. The older dense theoretical insertion-site track is
off by default; add `--show-theoretical-sites` if you want it too. Remove the
blue bars with `--no-read-histogram`.

For a single-contig reference, a CSV may omit `contig`; the GenBank accession
is used automatically. For multiple contigs, provide a contig column or map an
incoming name explicitly:

```bash
phage-tnseq-viz plot reference.gbk --final-dataset counts.csv \
  --contig-alias tpp_contig_1=NC_012345.1
```

### Final-site CSV schema

The canonical required columns are:

```csv
contig,position,read_count
NC_012345.1,287,14
NC_012345.1,415,0
```

- `position` is 1-based genomic insertion coordinate.
- `read_count` is the final count after any selected cutoff/averaging.
- Common imports are accepted: `accession`/`chrom` for `contig`,
  `TA_site`/`insertion_site` for `position`, and `count`/`mean_count` for
  `read_count`.
- Optional `raw_read_count`, `read_count_sd`/`sd_count`, and `n_subsamples`
  are preserved in the normalised output.
- Duplicate rows at the same coordinate are merged by summing counts.

For saturation-based essentiality, the table must represent **every candidate
insertion site, including zero-count sites**. By default the direct path fills
missing sites using the selected motif (`mariner`/`himar1` = `TA`). Set
`--candidate-model observed` only when you intentionally want to retain an
observed-sites-only table; it is not suitable for the built-in saturation
classifier. Available models are:

| Model | Meaning |
|---|---|
| `auto` | Motif candidates when a motif is known; all bases otherwise. |
| `motif` | All positions matching `--transposon` IUPAC motif. |
| `all-bases` | Every genomic base; a generic near-random-transposon heuristic. |
| `observed` | Do not create unobserved zero-count sites. |

The supplied essentiality thresholds were written for TA-site data. Applying
them to `all-bases`/near-random insertion data is an explicit adaptation, not
a biologically validated replacement for a transposon-specific method.

For example, a final table from a Tn5-like/near-random library should make that
choice explicit rather than silently inheriting the default mariner/TA model:

```bash
phage-tnseq-viz plot reference.gbk --final-dataset tn5_final.csv \
  --transposon tn5 --candidate-model all-bases
```

## Process raw Illumina reads

Run all selected stages non-interactively, which is best for reproducibility:

```bash
phage-tnseq-viz process phold_annotated_phage.gbk \
  --reads1 library_R1.fastq.gz --reads2 library_R2.fastq.gz \
  --output-dir run_phage_A \
  --quality-phred 25 --min-read-length 40 \
  --tpp-mode himar1 \
  --tpp-primer ACTTATCAGCCAACCTGTTA \
  --tpp-mismatches 1 \
  --min-mapped-reads 2
```

For a dual-IR cargo construct, replace `--tpp-primer` with the actual terminal
sequence observed at the transposon/genome junction. TPP's `Sassetti`/Himar1
and `Tn5` modes are exposed as `--tpp-mode himar1` and `--tpp-mode tn5`.
Check a known junction against the reference before interpreting biological
calls, especially when the library architecture differs from a standard TPP
protocol.

The command creates:

```text
run_phage_A/
  reference.fasta                 # GenBank converted for TPP/BWA
  01_fastqc/                       # if enabled
  02_fastp/                        # if enabled, including HTML/JSON report
  03_tpp/                          # TPP WIG / diagnostic files
  03_subsamples/                   # if depth matching is enabled
  pipeline.log
  processing_manifest.json         # exact commands, paths, parameters
  final_sites.csv
  final_gene_essentiality.csv
  tnseq_map.png
```

### Stage controls and interactive setup

Every external stage can be skipped. For example:

```bash
# Keep raw input unchanged but still run TPP+BWA
phage-tnseq-viz process reference.gbk --reads1 R1.fq.gz \
  --skip-fastqc --skip-fastp

# Only run the optional QC/filtering work; no final table/map is made
phage-tnseq-viz process reference.gbk --reads1 R1.fq.gz \
  --skip-tpp
```

To be prompted for the major stage choices, quality cutoff, length cutoff,
TPP mismatch allowance, count cutoff, and depth-matching settings, add
`--interactive`:

```bash
phage-tnseq-viz process reference.gbk --reads1 R1.fq.gz --reads2 R2.fq.gz --interactive
```

All answers have CLI equivalents and are written to `processing_manifest.json`;
interactive mode is a convenience, not a hidden configuration format.

### Read-count threshold and cross-phage depth matching

`--min-mapped-reads N` applies the requested strict rule: sites with fewer than
`N` mapped reads become zero; a count equal to `N` stays non-zero.

For fairer depth-matched comparisons between phage libraries, use the same
target read depth and number of replicates for each run:

```bash
phage-tnseq-viz process phage_A.gbk --reads1 A_R1.fq.gz --reads2 A_R2.fq.gz \
  --subsample-depth 1000000 --subsample-replicates 5 --subsample-seed 101

phage-tnseq-viz process phage_B.gbk --reads1 B_R1.fq.gz --reads2 B_R2.fq.gz \
  --subsample-depth 1000000 --subsample-replicates 5 --subsample-seed 101
```

SeqKit `sample2 -n -2` is called with the same seed for paired R1/R2 files,
then each subsample is independently processed by TPP. Counts are thresholded
per subsample and the final CSV reports their mean and population standard
deviation. The maps remain separate because phage genome coordinates are not
assumed to be alignable.

## Essentiality classification and colours

The built-in classifier is a faithful, dependency-free port of the supplied
`essentiality_classification.R` logic. Intergenic sites are retained in the
site CSV but excluded from gene calls; an overlapping site contributes to each
overlapping CDS.

| Rule | Built-in call |
|---|---|
| ≤5 candidate sites in a gene | No R call (`Insufficient sites` on map) |
| Saturation `< 0.20` | Essential |
| Saturation `> 0.80` | Non-essential |
| Otherwise | Strand-aware 3′-bias and first/second-half domain checks → Essential or Intermediate |
| Saturation 0.20–0.40; >50% positive sites below pooled positive-count median | Strong fitness defect (unless Essential) |
| Saturation 0.70–0.80; >50% positive sites above that median | Reduced fitness (unless Essential) |

Saturation is `positive candidate sites / all candidate sites`. Arrow fills use
red (Essential), orange (Strong fitness defect), pale yellow (Intermediate),
light blue (Reduced fitness), green (Non-essential), and grey for
insufficient/ambiguous calls. Arrow borders continue to use the phold/pharokka
PHROG functional-category palette, so functional annotation and fitness call
remain independently readable.

### Write a custom classifier

Copy [examples/custom_classifier_template.py](examples/custom_classifier_template.py)
and change its function:

```python
def classify_gene(gene_id, site_rows):
    # site_rows is a tuple of GeneSiteRow objects:
    # contig, gene_id, strand, position, read_count
    return "My custom call"  # or None
```

Then pass it to either path:

```bash
phage-tnseq-viz plot reference.gbk --final-dataset final.csv \
  --classifier my_classifier.py
```

The classifier file is ordinary local Python and executes with your user
permissions. Only use a file you trust. Unrecognised custom labels are retained
verbatim in the CSV and shown in neutral grey in the built-in plot legend.

## Legacy pre-visualisation mode

All pre-v0.2 commands still work:

```bash
phage-tnseq-viz genome.gbk -t mariner -o phage_map.png
phage-tnseq-viz genome.gbk --gc-content --gc-skew --trna
phage-tnseq-viz genome.gbk --wrap-kb 10 --paper a4 --fit-page
```

Use `--no-wrap`, `--rows`, `--small-title`, `--transparent`, `--dpi`, and the
other existing map-style options as before. Run `phage-tnseq-viz --help`,
`phage-tnseq-viz plot --help`, or `phage-tnseq-viz process --help` for every
parameter.

## Verification

The test suite uses synthetic data and command-plan tests; it never requires
FastQC, fastp, TRANSIT, BWA, or SeqKit binaries:

```bash
pytest
```

## License

See [LICENSE](LICENSE).
