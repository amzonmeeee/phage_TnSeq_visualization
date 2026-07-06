# Phage Tn-Seq Visualization Tool

**Pre-visualization** for phage transposon-sequencing (Tn-Seq / HIDDEN-seq) experiments.

Given a phage genome, this tool draws a **linear genome map** with gene arrows and
the positions where a chosen transposon can insert — *before* any sequencing data is
available. It is designed to fit naturally next to the
[pharokka](https://github.com/gbouras13/pharokka) →
[phold](https://github.com/gbouras13/phold) →
[phynteny](https://github.com/susiegriggo/Phynteny_transformer) annotation pipeline,
and reuses their PHROG functional-category colour scheme.

> This is v0.1 (pre-visualization). Sequencing-data / essentiality overlays are planned
> for a later version — the gene-arrow **fill** colour is reserved for that.

## What it draws

A horizontal genome line (length-scaled, in kb) with, per genome/contig:

* **Gene / ORF arrows** — arrow direction shows strand.
  * **Fill = essentiality** → drawn **black** for now (no sequencing data yet).
  * **Border colour = PHROG functional category** (phold/pharokka colour scheme, with legend).
  * De-novo ORFs (from an unannotated genome) get a **grey dashed border**.
* **Insertion-site track** — short **red** ticks at every position the transposon can insert.
* Optional **insertion-density** heat track (sites/kb, colour-mapped) to spot hot/cold zones.
* Optional **GC content**, **GC skew**, and **tRNA/CRISPR** tracks.
* Genome **name / accession** label (auto-detected or user-supplied).

## Inputs

Two flavours of GenBank (`.gbk` / `.gb`) are accepted:

1. **Annotated** — output of pharokka / phold / phynteny(_transformer). CDS `function`
   qualifiers (PHROG categories) drive the border colours.
2. **Plain / unannotated** — e.g. a GenBank download. ORFs are called de-novo with
   [pyrodigal-gv](https://github.com/althonos/pyrodigal-gv) (the same viral gene caller
   pharokka uses internally).

## Built on existing tools (no reinvented wheels)

| Job | Library |
|-----|---------|
| GenBank parsing | [Biopython](https://biopython.org) |
| De-novo ORF calling | [pyrodigal-gv](https://github.com/althonos/pyrodigal-gv) |
| Linear genome plotting | [pyGenomeViz](https://github.com/moshi4/pyGenomeViz) (same author as pyCirclize, which phold/pharokka use) |
| Insertion-site scanning | Biopython + IUPAC motif regex |

## Installation

Pure `pip`, no conda/homebrew required (works on system Python ≥ 3.9):

```bash
git clone <this-repo>
cd phage_TnSeq_visualization
pip install .
```

For development:

```bash
pip install -e ".[dev]"
python examples/make_example.py   # generate synthetic test genomes
pytest
```

## Usage

```bash
phage-tnseq-viz INPUT.gbk -t mariner -o phage_map.png
```

### Transposon insertion preference

`-t / --transposon` accepts either a **preset** or a **custom IUPAC motif**:

| Preset | Motif |
|--------|-------|
| `mariner`, `himar1` | `TA` |
| `tn5`, `tn10`, `tn7`, `mu` | ≈ random / site-specific (ticks skipped) |

```bash
phage-tnseq-viz genome.gbk -t mariner        # TA sites
phage-tnseq-viz genome.gbk -t NTAN           # custom IUPAC motif
phage-tnseq-viz genome.gbk -t TA --single-strand
```

### Common options

```bash
# Extra tracks
phage-tnseq-viz genome.gbk --insertion-density        # sites/kb heat track
phage-tnseq-viz genome.gbk --gc-content --gc-skew --trna

# Naming & legend
phage-tnseq-viz genome.gbk --name "Phage vB_EcoM_XYZ"
phage-tnseq-viz genome.gbk --no-legend

# Output size / paper / format / transparency
phage-tnseq-viz genome.gbk -o map.svg                     # SVG (vector)
phage-tnseq-viz genome.gbk --paper a4 --fit-page          # fit A4 landscape
phage-tnseq-viz genome.gbk --paper a3 --portrait
phage-tnseq-viz genome.gbk --width 20 --track-height 2    # custom figure size
phage-tnseq-viz genome.gbk --transparent --dpi 600
```

Run `phage-tnseq-viz --help` for the full list.

### Output

PNG (raster, `--dpi` configurable) or SVG (vector) — chosen by the output file
extension. `--transparent` gives a transparent background for either.

## Roadmap

* Overlay real Tn-Seq read counts → colour arrow **fill** by gene essentiality.
* Per-insertion-site read-count heat track.
* Multi-genome comparison / alignment view.

## License

See [LICENSE](LICENSE).
