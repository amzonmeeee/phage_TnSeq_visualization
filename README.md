# Phage Tn-Seq Visualization Tool

**Pre-visualization** for phage transposon-sequencing (Tn-Seq / HIDEN-seq) experiments.

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

A horizontal genome line (length-scaled, in kb) that **wraps onto multiple rows**
(≈ every 20 kb by default, evenly balanced) so both short and long genomes stay
paper-friendly — tune with `--wrap-kb`, force a row count with `--rows`, or draw a
single line with `--no-wrap`. Per genome/contig it shows:

* **Gene / ORF arrows** — arrow direction shows strand.
  * **Fill = essentiality** → drawn **black** for now (no sequencing data yet).
  * **Border colour = PHROG functional category** (phold/pharokka colour scheme, with legend).
  * De-novo ORFs (from an unannotated genome) get a **grey dashed border**.
* **Insertion-site track** — short **red** ticks at every position the transposon can insert.
* **Insertion-density** heat track (sites/kb, colour-mapped) to spot hot/cold zones; on by default.
* Optional **GC content**, **GC skew**, and **tRNA/CRISPR** tracks.
* Genome **name / accession** label (auto-detected or user-supplied).

## Inputs

Two flavours of GenBank (`.gbk` / `.gb`) are accepted:

1. **Annotated** — output of pharokka / phold / phynteny(_transformer). CDS `function`
   qualifiers (PHROG categories) drive the border colours.
2. **Plain / unannotated** — e.g. a GenBank download. ORFs are called de-novo with
   [pyrodigal-gv](https://github.com/althonos/pyrodigal-gv) (the same viral gene caller
   pharokka uses internally).

## Built on existing tools

| Job | Library |
|-----|---------|
| GenBank parsing | [Biopython](https://biopython.org) |
| De-novo ORF calling | [pyrodigal-gv](https://github.com/althonos/pyrodigal-gv) |
| Linear genome plotting | [pyGenomeViz](https://github.com/moshi4/pyGenomeViz) (same author as pyCirclize, which phold/pharokka use) |
| Insertion-site scanning | Biopython + IUPAC motif regex |

## Installation (pending review)

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

E.g.
```bash
phage-tnseq-viz genome.gbk -t mariner        # TA sites
phage-tnseq-viz genome.gbk -t NTAN           # custom IUPAC motif
phage-tnseq-viz genome.gbk -t TA --single-strand
```

### Common options

```bash
# Extra tracks
phage-tnseq-viz genome.gbk --no-insertion-density     # drop the sites/kb heat track
phage-tnseq-viz genome.gbk --gc-content --gc-skew --trna

# Naming & legend
phage-tnseq-viz genome.gbk --name "Phage vB_EcoM_XYZ"
phage-tnseq-viz genome.gbk --small-title          # name as a small side label, not a heading
phage-tnseq-viz genome.gbk --no-legend

# Line wrapping (multi-row layout for papers)
phage-tnseq-viz genome.gbk                       # auto-wrap ~every 20 kb
phage-tnseq-viz genome.gbk --wrap-kb 10          # new row about every 10 kb
phage-tnseq-viz genome.gbk --rows 3              # force exactly 3 rows
phage-tnseq-viz genome.gbk --no-wrap             # single long line

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

## Citation
blahblah

## License

See [LICENSE](LICENSE).
