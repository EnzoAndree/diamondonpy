# DIAMOND Wrapper Real-World Usage Scenarios

This directory contains examples of how to use the `diamondonpy` Python wrapper for DIAMOND in real-world bioinformatics applications.

## Overview

The `real_world_scenarios.py` script demonstrates five common bioinformatics workflows using the DIAMOND aligner:

1. **Genome Annotation Pipeline** - Using BLASTX for annotating newly sequenced genomes
2. **Protein Family Analysis** - Clustering proteins into families based on sequence similarity
3. **Ortholog Identification** - Finding orthologous genes between species using bidirectional best hits
4. **Metagenomic Analysis** - Analyzing the taxonomic composition of metagenomic samples
5. **Protein Domain Analysis** - Identifying conserved domains in protein sequences

## Requirements

- DIAMOND (command-line tool should be installed and in PATH)
- diamondonpy
- pandas
- matplotlib (for visualizations)
- seaborn (optional, for enhanced visualizations)

## Usage

You can run individual scenarios:

```bash
# Run a specific scenario (1-5)
python real_world_scenarios.py --scenario 1

# Run all scenarios
python real_world_scenarios.py --all
```

## Detailed Scenario Descriptions

### Scenario 1: Genome Annotation Pipeline

This scenario demonstrates how to annotate coding sequences from a newly sequenced genome against a reference protein database.

**Steps:**
1. Download or create a reference protein database
2. Create a DIAMOND database 
3. Use BLASTX to translate nucleotide sequences and align against the reference
4. Filter and analyze the results
5. Generate an annotation report

**Key Features:**
- Handles both nucleotide and protein sequences
- Performs six-frame translation (via DIAMOND blastx)
- Outputs annotation results as a structured DataFrame
- Demonstrates filtering by bit score, E-value, and sequence identity

### Scenario 2: Protein Family Analysis through Clustering

This scenario shows how to group proteins into families or clusters based on sequence similarity.

**Steps:**
1. Create or use a dataset of related proteins
2. Run DIAMOND clustering to group similar sequences
3. Analyze the resulting clusters
4. Visualize cluster sizes and relationships

**Key Features:**
- Uses DIAMOND's efficient clustering algorithms
- Groups proteins based on sequence similarity
- Identifies representative sequences for each cluster
- Visualizes cluster sizes for easy interpretation

### Scenario 3: Ortholog Identification

This scenario demonstrates how to identify orthologous genes between two species using the bidirectional best hit approach.

**Steps:**
1. Create protein datasets for two species
2. Create DIAMOND databases for both species
3. Perform bidirectional best hit analysis
4. Identify and analyze ortholog pairs

**Key Features:**
- Performs reciprocal best-hit searches
- Filters results by E-value and sequence identity
- Identifies ortholog pairs between species
- Reports sequence similarity metrics for each pair

### Scenario 4: Metagenomic Analysis

This scenario shows how to analyze metagenomic data to identify the taxonomic composition of a microbial community.

**Steps:**
1. Process metagenomic contigs or reads
2. Translate and align against a reference protein database
3. Analyze taxonomic assignments
4. Generate visualizations of community composition

**Key Features:**
- Processes assembled contigs or individual reads
- Matches sequences to taxonomic information
- Analyzes community composition
- Visualizes taxonomic distribution

### Scenario 5: Protein Domain Analysis

This scenario demonstrates how to identify conserved domains in protein sequences by alignment to known domain databases.

**Steps:**
1. Prepare protein sequences of interest
2. Align against a domain database
3. Identify and map domains to protein sequences
4. Visualize domain architecture

**Key Features:**
- Identifies conserved domains in proteins
- Maps domains to specific regions in proteins
- Creates simple domain architecture visualizations
- Detects multiple domains within a single protein

## Data Storage

The script creates a `data` directory to store:
- Input FASTA files
- DIAMOND databases
- Analysis results
- Generated visualizations

## Extending the Examples

These examples provide a foundation for building more complex bioinformatics pipelines. You can extend them by:

1. Using larger, real-world datasets
2. Integrating with other bioinformatics tools
3. Adding more sophisticated analysis methods
4. Creating more detailed visualizations
5. Implementing parallel processing for large-scale analyses

## References

- DIAMOND: [https://github.com/bbuchfink/diamond](https://github.com/bbuchfink/diamond)
- diamondonpy: [https://github.com/yourusername/diamondonpy](https://github.com/yourusername/diamondonpy) 