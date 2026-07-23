# Phage Tn-Seq Visualization Tool

`phage-tnseq-viz` turns an annotated phage GenBank reference and Tn-Seq data
into a linear genome map, a per-insertion-site CSV, and gene-level
essentiality calls. It can either orchestrate optional raw-read processing or
start directly from a final count table produced elsewhere.

The tool is designed for transposon-enriched Illumina libraries after sheared
genomes, adapter ligation, transposon-specific PCR, and junction sequencing.
It does not attempt to infer a non-standard IR–cargo–IR primer sequence which is is an experimental input and must be supplied to TPP when needed.

## What is new in v0.2

- Optional FastQC → fastp → TRANSIT TPP + BWA → SeqKit preprocessing path.
- A direct `--final-dataset` path that requires no external bioinformatics
  programs.
- Blue per-site read-count bars above gene arrows.
- Essentiality-coloured gene fills, while PHROG functional categories remain
  arrow borders.
- Normalised site-level and gene-level CSV output, processing logs, and a
  reproducibility manifest.
- A dependency-free implementation of the supplied Harms-lab R classification<sup>1</sup>
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
git clone https://github.com/amzonmeeee/phage_TnSeq_visualization
cd phage_TnSeq_visualization
pip install .
```

For development:

```bash
pip install -e ".[dev]"
pytest
```

The Python package includes Biopython, pyGenomeViz, matplotlib, and
pyrodigal-gv. Interactive `.html` maps additionally need Plotly, installed via
the optional `[interactive]` extra (`pip install ".[interactive]"`). The
raw-read route deliberately does **not** install external programs for you. Install the tools you plan to run (normally in a dedicated
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

This writes `example_map.png`, `example_map_sites.csv`,
`example_map_gene_essentiality.csv`, and `example_map_qc.csv`.

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
and essentiality gene arrows. The older dense theoretical insertion-site track is
off by default; add `--show-theoretical-sites` if you want it too. Remove the
blue bars with `--no-read-histogram`.

Real Tn-Seq counts are heavy-tailed, so a single hypersaturated site can squash
every other blue bar. Cap the read-count scale at a percentile of each contig's
positive counts with `--read-histogram-cap 95`; taller bars clip to full height
and the scale then reads `≥ N`.

### Interactive HTML map

The output format follows the `--output` extension: `.png`/`.svg` are the static
maps, and **`.html` renders an interactive, zoomable map** (pan/zoom, and hover
to read the exact per-site count, gene product, PHROG category, and essentiality
call). The map is built with [Plotly](https://plotly.com/python/), and its
`plotly.js` library is embedded directly in the file — the HTML is therefore
fully self-contained and opens offline in any browser with no network access
(and so each `.html` map is a few MB).

```bash
phage-tnseq-viz plot phage.gbk --final-dataset counts.csv -o phage_map.html
```

Interactive output needs Plotly, kept as an optional extra so the static path
stays lightweight:

```bash
pip install ".[interactive]"   # or: pip install plotly
```

### Progress verbosity

Every command prints stage progress by default. Use `-q`/`--quiet` to print only
warnings and errors (handy in scripts and batch runs), or `-v`/`--verbose` for
extra detail. Progress goes to stdout, warnings/errors to stderr.

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
  final_qc.csv                    # library QC metrics
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

All answers have CLI equivalents and are written to `processing_manifest.json`.
Interactive mode is a convenience, not a hidden configuration format.

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

### TTR normalization (an alternative to subsampling)

Subsampling makes libraries comparable by throwing reads away until they share a
depth. `--normalize ttr` reaches the same goal without discarding anything: it
rescales each contig so its typical non-zero count meets a fixed target (default
100), so a library run entirely on its own still lands on the same scale as any
other.

```bash
phage-tnseq-viz plot phage_A.gbk --final-dataset A_final.csv -o A.png --normalize ttr
phage-tnseq-viz plot phage_B.gbk --final-dataset B_final.csv -o B.png --normalize ttr
```

The factor is `target / (density × trimmed_mean_of_hit_sites)`, dropping the top
and bottom 5% of non-zero counts before averaging so a few hypersaturated sites
cannot set the scale for the whole library. This is TRANSIT's default (TTR =
Trimmed Total Reads), available on both the `plot` and `process` paths and
computed per contig, since phage contigs are not sequenced to comparable depths.
Change the target with `--norm-target`.

Two things are worth being clear about:

- **It never changes a within-library call.** Scaling a whole contig by one
  number leaves saturation untouched and moves every count together, so the
  essentiality calls are identical with or without it. Its effect is entirely on
  *cross*-library comparability and on the read-count scale of the map.
- **QC still reports raw counts.** Normalization pins the non-zero mean to the
  target, which would otherwise hide the low-`NZmean` signature of an
  under-sequenced library, so the QC table is always computed before scaling.
  The normalized value goes in the `read_count` column and the original is kept
  in `raw_read_count`.

TTR corrects for depth only. It does not fix a position-dependent coverage
gradient or a genuinely different insertion-site distribution between libraries.

## Library quality control

Every essentiality call rests on assumptions about the library: that enough
candidate sites were hit to tell an essential gene from an unlucky one, and that
the counts are not dominated by a few runaway sites. Both paths therefore print
a QC table and write it to `<stem>_qc.csv`, so a result can be judged before it
is believed.

```text
     contig  sites  hit  density   NZmean  NZmedian      max  skew   PTI
  NC_TEST01   1293  188    0.145  21280.1       4.0  4000000  13.6  0.00
warning: library QC (NC_TEST01): saturation 0.15 is below 0.30
warning: library QC (NC_TEST01): top site has 4000000 reads, a likely outlier
```

That example shows why `NZmedian` sits next to `NZmean`: a mean of 21,280 against
a median of 4 is the signature of a single site swallowing the library, which is
also what `--read-histogram-cap` exists to make the map readable despite.

| Metric | Meaning |
|---|---|
| `density` | Saturation: the fraction of candidate sites that were hit. |
| `mean_count` | Mean count over **all** sites, zeros included. |
| `nz_mean` / `nz_median` | Mean and median over hit sites only; a wide gap means outliers. |
| `max_count` | Largest single-site count. |
| `skewness` / `kurtosis` | Third and fourth moments over hit sites (kurtosis is excess, so 0 is normal). |
| `pickands_tail_index` | Tail-weight estimate; higher means a few sites carry disproportionate reads. |

The metric set and the warning thresholds follow TRANSIT's `tnseq_stats`, which
is the established vocabulary for Tn-Seq QC. Two things differ, both because
this tool targets phage:

- Metrics are reported **per contig**, matching how saturation and read
  thresholds are already computed elsewhere here.
- The Pickands tail index is **adapted to small genomes**. It compares order
  statistics at ranks M, 2M and 4M; TRANSIT scans a fixed M = 10..99, which
  silently requires at least 397 sites and raises `IndexError` below that — and
  a phage genome routinely has fewer. Here the scan is capped by the data
  instead, so it agrees with TRANSIT exactly whenever TRANSIT's range applies
  and reports nothing rather than crashing when it does not. Ranks that tie,
  which happens once a sparse library's tail runs into the zero counts, are
  skipped rather than poisoning the median.

The thresholds themselves were derived from bacterial libraries. On a phage
genome, with far fewer candidate sites, treat a warning as a prompt to look
rather than a pass/fail verdict.

## Essentiality classification and colours

The built-in classifier is a faithful, dependency-free port of the supplied
`essentiality_classification.R` logic<sup>1</sup>. Intergenic sites are retained in the
site CSV but excluded from gene calls; an overlapping site contributes to each
overlapping CDS.

| Rule | Built-in call |
|---|---|
| ≤5 candidate sites in a gene | No R call (`Insufficient sites` on map, unless the binomial supplement applies) |
| Saturation `< 0.20` | Essential |
| Saturation `> 0.80` | Non-essential |
| Otherwise | Strand-aware 3′-bias and first/second-half domain checks → Essential or Intermediate |
| Saturation 0.20–0.40; >50% positive sites below pooled positive-count median | Strong fitness defect (unless Essential) |
| Saturation 0.70–0.80; >50% positive sites above that median | Reduced fitness (unless Essential) |

Saturation is `positive candidate sites / all candidate sites`. Arrow fills use
maroon (Essential), bright red (Essential — binomial), red (Strong fitness
defect), orange (Intermediate), amber (Reduced fitness), cream (Non-essential),
and grey for insufficient/ambiguous calls. Arrow borders continue to use the
phold/pharokka PHROG functional-category palette, so functional annotation and
fitness call remain independently readable.

### The binomial supplement for short genes

The R rules refuse to call any gene with five or fewer candidate sites. On a
phage genome that silences a real share of the annotation, because genes are
short and TA sites are sparse — especially in GC-rich phages, where a 200 bp ORF
may hold only three or four TA sites. Yet a gene with *no* insertions at all is
the clearest essentiality signal there is, and whether that is surprising depends
on how saturated the library is rather than on a fixed site count.

So a gene the R rules leave uncalled is re-examined: with a library-wide
non-insertion probability `phi = 1 - saturation`, a gene of `n` candidate sites
is empty by chance with probability `phi**n`. When that is below 0.05 — that is,
when `n > log(0.05)/log(phi)` — the gene is called **`Essential (binomial)`**.
The threshold is reported per contig in the `binomial_min_sites` column of the
gene CSV, so every such call can be checked by hand.

This is a supplement, not a replacement: it only ever fills a result the R rules
left empty, and can never change a call they made. It follows the binomial
addition Choudhery et al. made to TRANSIT's Gumbel method<sup>2</sup>, where these
genes are reported as `EB` alongside Gumbel's `E`. It scales with library
quality by construction — a well-saturated library recovers genes down to two or
three sites, while a sparse one recovers none, which is the conservative and
correct outcome.

Turn it off with `--no-binomial-short-genes` to reproduce the R script exactly:

```bash
phage-tnseq-viz plot reference.gbk --final-dataset final.csv \
  --no-binomial-short-genes
```

### Gap analysis: a second opinion the saturation rules cannot give

Saturation counts *how many* sites in a gene were hit, never *where* they were.
Two genes with identical saturation therefore score identically, even when one
has its empty sites scattered evenly and the other has them piled into a single
unbroken block. The second is the interesting one, because an uninterrupted
stretch of dead sites is what a domain that cannot tolerate insertion looks
like.

Every gene is therefore also scored on its **longest run of consecutive empty
sites**. Under a model where each candidate site is hit independently with the
library-wide probability, that longest run follows a Gumbel (extreme-value)
distribution, so the observed run gets a p-value. P-values are corrected across
all genes with Benjamini-Hochberg, and the gene CSV reports `max_gap`,
`expected_max_gap`, `gap_pvalue` and `gap_qvalue`.

This is the method of Griffin et al. (2011)<sup>3</sup>, with the refinement
TRANSIT later adopted: the distribution's location is fixed by matching moments
against the expected longest run, which is computed *exactly* (via Boyd's
recurrence) for genes under 20 sites. That exactness is not a nicety here —
phage genes are mostly that small, and the usual asymptotic formula is poor in
that range.

**The gap analysis never changes an essentiality call.** It is independent
evidence, and its value is in the disagreements. When a gene's dead run is
significant but the rules did *not* call it essential, the `gap_flag` column
reads `domain-essential?`. That is the multi-domain case the R rules
structurally cannot catch: they split a gene at its midpoint and ask whether one
half is >70% empty and the other <30%, so a dead domain lying across that split
defeats the test and the gene escapes as `Intermediate`. Treat the flag as a
prompt to look at the gene, not as a call.

Skip the whole analysis with `--no-gap-analysis`; the columns are then left
blank and every other output is unchanged.

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

## References

1. Humolli D, Ransome J, Piel D, Veening JW, Harms A. Systematic mapping of bacteriophage gene essentiality with HIDEN‑SEQ. _bioRxiv (Cold Spring Harbor Laboratory)_. Published online November 20, 2025. doi:10.1101/2025.11.20.689424
2. Choudhery S, Brown AJ, Akusobi C, Rubin EJ, Sassetti CM, Ioerger TR. Modeling site-specific nucleotide biases affecting Himar1 transposon insertion frequencies in TnSeq data sets. _mSystems_. 2021;6(5):e00876-21. doi:10.1128/mSystems.00876-21
3. Griffin JE, Gawronski JD, DeJesus MA, Ioerger TR, Akerley BJ, Sassetti CM. High-resolution phenotypic profiling defines genes essential for mycobacterial survival and cholesterol catabolism. _PLoS Pathogens_. 2011;7(9):e1002251. doi:10.1371/journal.ppat.1002251
4. Schilling MF. The longest run of heads. _College Mathematics Journal_. 1990;21(3):196-207. doi:10.1080/07468342.1990.11973306

## Citation

If you want to cite this tool, please site blahblahblah.

## License

See [LICENSE](LICENSE).
