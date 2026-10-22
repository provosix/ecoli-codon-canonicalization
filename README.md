# Synonymous Codon Canonicalization in *E. coli* K-12 MG1655

An independent bioinformatics project that compares three strategies for replacing synonymous codons in annotated coding DNA sequences (CDSs) with one selected codon per amino acid.

## Project question

How do different synonymous-codon selection strategies change coding-sequence composition in *E. coli* K-12 MG1655?

The workflow compares three models:

- **Most frequent:** choose the most frequently observed synonymous codon for each amino acid.
- **Random:** choose one synonymous codon per amino acid using a fixed random seed (`42` by default).
- **Least frequent:** choose the least frequently observed synonymous codon for each amino acid.

For each model, the script measures codon and nucleotide changes, GC content, and changes by codon position. It also checks that each redesigned CDS translates to the same protein sequence as its original annotated CDS.

## Data and methods

- **Organism/strain:** *Escherichia coli* K-12 MG1655
- **GenBank accession:** [`NC_000913.3`](https://www.ncbi.nlm.nih.gov/nuccore/NC_000913.3)
- **Input file:** `data/ecoli.gb`
- **Translation table:** NCBI bacterial/plastid code 11
- **Random seed:** `42`

The script excludes pseudogene annotations, CDSs with selenocysteine translation exceptions, CDSs annotated with ribosomal slippage, and CDSs that fail sequence validation. The terminal stop codon is preserved in FASTA output but excluded from sense-codon composition statistics.

The first codon of each redesigned CDS is set to `ATG` to preserve the start-methionine convention used by this analysis. Alternative bacterial initiation codons require special handling because their translation depends on their position at the start of a CDS.

## Results from the included run

| Model | Sense codons changed | Sense nucleotides changed | GC after canonicalization |
|---|---:|---:|---:|
| Most frequent | 48.29% | 18.57% | 57.56% |
| Random (seed 42) | 65.92% | 23.65% | 46.08% |
| Least frequent | 77.03% | 29.35% | 45.36% |

These are computational sequence-composition results. They do **not** demonstrate that redesigned sequences preserve expression levels, RNA structure, translation dynamics, or organismal fitness in living cells.

## Requirements

- Python 3.10 or newer
- Biopython
- Matplotlib
- openpyxl

Install dependencies from the project root:

```bash
python -m pip install -r requirements.txt
```

## Run the analysis

From the project root:

```bash
python analyze_ecoli.py --input data/ecoli.gb --output-dir results_generated
```

To use a different random seed:

```bash
python analyze_ecoli.py --input data/ecoli.gb --output-dir results_generated --seed 123
```

## Output files

The output directory contains:

- `main_results.xlsx` — main comparison, dataset details, codon usage, and model mappings.
- `supplementary_tables.xlsx` — codon-level and amino-acid-level tables.
- `results_summary.json` — machine-readable summary of results.
- `model_comparison.png` — comparison plot generated from computed values.
- `canonicalized_cds_most_frequent.fasta` — redesigned CDS sequences for the most-frequent model.
- `canonicalized_cds_random.fasta` — redesigned CDS sequences for the random model.
- `canonicalized_cds_least_frequent.fasta` — redesigned CDS sequences for the least-frequent model.

## Reproducibility and limitations

- Results are calculated from the GenBank input; the plot and summary tables are generated from computed values.
- The random model is reproducible for the same input, software dependencies, and seed.
- The script validates protein-sequence preservation for retained CDSs and stops if a mismatch is detected.
- This is a sequence-level computational project, not experimental validation of redesigned genes.

## Key Findings

The analysis compared three synonymous codon canonicalization strategies in *E. coli* K-12 MG1655.

| Strategy | Sense codons changed | Nucleotides changed | Final GC content |
|---|---:|---:|---:|
| Most frequent codon | 48.29% | 18.57% | 57.56% |
| Random codon | 65.92% | 23.65% | 46.08% |
| Least frequent codon | 77.03% | 29.35% | 45.36% |

### Main observations

- The most-frequent-codon strategy produced the lowest proportion of changed sense codons and nucleotides among the three strategies.
- The least-frequent-codon strategy produced the highest proportion of changes.
- GC content varied across strategies, demonstrating that synonymous codon selection can substantially alter nucleotide composition.

These findings describe computational sequence comparisons. They do not establish equivalent gene expression, translation efficiency, or biological fitness in living cells.
