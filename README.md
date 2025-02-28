# 💎 Diamononpy 🐍

<img src="img/diamondonpy.png" align="left" width="190" style="padding-right: 20px">

A Python wrapper for the ultra-fast [DIAMOND](https://github.com/bbuchfink/diamond/) sequence alignment tool. This package provides a clean, Pythonic API for DIAMOND's powerful sequence search capabilities with seamless pandas integration for efficient bioinformatics data analysis and processing. Perfect for researchers and bioinformaticians working with large genomic datasets who need both speed and ease of use.

### ✨ Features
- 🚀 Full support for all DIAMOND V commands
- 📊 Results returned as pandas DataFrames for easy analysis
- 🗑️ Automatic temporary file management
- 🔍 Type hints for better IDE support
- 🧪 Comprehensive test suite
- 📦 Minimal dependencies (only pandas and numpy)

<br clear="left">

## 📥 Installation

First, ensure you have DIAMOND installed and accessible in your PATH. Then install this package:

```bash
pip install diamononpy
```

For installation directly from the GitHub repository, use the following command:

```bash
pip install git+https://github.com/EnzoAndree/diamondonpy.git
```

For development installation (including test dependencies):

```bash
pip install -e ".[dev]"
```

## 🚀 Usage

### Basic Usage

```python
from diamononpy import Diamond

# Initialize the wrapper
diamond = Diamond()

# Create a database
diamond.makedb(
    db="mydb.dmnd",
    input_file="sequences.fasta",
    threads=4
)

# Run BLASTP search - results as DataFrame
results_df = diamond.blastp(
    db="mydb.dmnd",
    query="query.fasta",
    evalue=1e-10,
    threads=4
)

# Access results using pandas
print(results_df.head())
print(f"Found {len(results_df)} hits")
print(f"Average identity: {results_df['pident'].mean():.2f}%")

# Filter results
significant_hits = results_df[
    (results_df['evalue'] < 1e-30) & 
    (results_df['pident'] > 90)
]
```

### 📊 Working with Results

All BLAST-like commands (blastp, blastx) return pandas DataFrames with the following columns:
- qseqid: Query sequence identifier
- sseqid: Subject sequence identifier
- pident: Percentage of identical matches
- length: Alignment length
- mismatch: Number of mismatches
- gapopen: Number of gap openings
- qstart: Start of alignment in query
- qend: End of alignment in query
- sstart: Start of alignment in subject
- send: End of alignment in subject
- evalue: Expect value
- bitscore: Bit score

```python
# BLASTP with output file
results_df = diamond.blastp(
    db="mydb.dmnd",
    query="query.fasta",
    out="results.txt",  # Optional: save to file
    evalue=1e-10
)

# Clustering with results as DataFrame
clusters_df = diamond.cluster(
    db="mydb.dmnd",
    approx_id=90.0
)
print(clusters_df.head())

# Bidirectional Best Hit analysis
bbh_df = diamond.bidirectional_best_hit(
    db1="db1.dmnd",
    db2="db2.dmnd",
    evalue=1e-10
)
print(bbh_df.head())
```

## 🛠️ Available Commands

All major DIAMOND commands are supported with enhanced result handling:

- `makedb`: Build DIAMOND database from a FASTA file
- `blastp`: Align protein sequences (returns DataFrame)
- `blastx`: Align DNA sequences (returns DataFrame)
- `view`: View DAA files (returns DataFrame for tabular output)
- `cluster`: Cluster sequences (returns DataFrame)
- `linclust`: Linear-time clustering (returns DataFrame)
- `getseq`: Retrieve sequences
- `dbinfo`: Database information
- `bidirectional_best_hit`: Perform bidirectional best hit analysis between two databases

## 🧠 Advanced Features

### Custom Output Formats

```python
# Custom BLAST output format
results_df = diamond.blastp(
    db="mydb.dmnd",
    query="query.fasta",
    outfmt="6 qseqid sseqid pident evalue bitscore qcovhsp"
)

# Non-tabular output
text_output = diamond.view(
    daa="alignment.daa",
    outfmt=0  # BLAST pairwise format
)
```

### 🧹 Temporary File Management

The package automatically manages temporary files:

```python
with Diamond() as diamond:
    results_df = diamond.blastp(
        db="mydb.dmnd",
        query="query.fasta"
    )
    # Temporary files are automatically cleaned up
```

### 🔄 Cluster Analysis

Analyze clustering results with built-in parser:

```python
# Perform clustering
clusters_df = diamond.cluster(
    db="mydb.dmnd",
    approx_id=90.0
)

# Data contains cluster IDs and members
print(f"Number of clusters: {clusters_df['cluster_id'].nunique()}")
```

## 📋 Requirements

- Python ≥ 3.6
- DIAMOND (must be installed separately and accessible in PATH)
- pandas ≥ 1.0.0
- numpy ≥ 1.18.0

## 🧪 Development

First, install development dependencies:

```bash
pip install -e ".[dev]"
```

To run tests:

```bash
# Run basic tests
pytest

# Run tests with coverage report
pytest --cov=diamononpy

# Run tests with detailed coverage report
pytest --cov=diamononpy --cov-report=term-missing

# Run tests verbosely
pytest -v

# Run a specific test file
pytest tests/test_diamond.py

# Run a specific test function
pytest tests/test_diamond.py::test_blastp
```

## 📚 References

This package is a wrapper for the [DIAMOND](https://github.com/bbuchfink/diamond/) bioinformatics tool.

When using DIAMOND in published research, please cite:

-   Buchfink B, Reuter K, Drost HG, \"Sensitive protein alignments at tree-of-life
    scale using DIAMOND\", *Nature Methods* **18**, 366–368 (2021).
    [doi:10.1038/s41592-021-01101-x](https://doi.org/10.1038/s41592-021-01101-x)

For sequence clustering:

-   Buchfink B, Ashkenazy H, Reuter K, Kennedy JA, Drost HG, \"Sensitive clustering
    of protein sequences at tree-of-life scale using DIAMOND DeepClust\", *bioRxiv*
    2023.01.24.525373; doi: [doi:10.1101/2023.01.24.525373](https://doi.org/10.1101/2023.01.24.525373)

Original publication to cite DIAMOND until v0.9.25:

-   Buchfink B, Xie C, Huson DH, \"Fast and sensitive protein alignment
    using DIAMOND\", *Nature Methods* **12**, 59-60 (2015).
    [doi:10.1038/nmeth.3176](https://doi.org/10.1038/nmeth.3176)

## 📄 License

This project is licensed under the GNU General Public License v3 (GPLv3) - see the [LICENSE](/LICENSE) file for details.