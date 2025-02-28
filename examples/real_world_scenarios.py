#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Real-world usage scenarios for the diamondonpy package.

This script demonstrates practical applications of diamondonpy in common bioinformatics workflows.
Each scenario is implemented as a function that can be run independently.

Requirements:
    - diamondonpy
    - pandas
    - matplotlib
    - seaborn (optional, for visualization)

Example usage:
    python real_world_scenarios.py --scenario 1

Author: [Your Name]
License: MIT
"""

import os
import sys
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from diamondonpy import Diamond
import tempfile
import urllib.request
import gzip
import shutil
from pathlib import Path

# Set up basic logging
import logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Define paths for data downloads
DATA_DIR = Path("data")
DATA_DIR.mkdir(exist_ok=True)

def download_file(url, output_path):
    """Download a file from a URL to the specified path."""
    logger.info(f"Downloading {url} to {output_path}")
    with urllib.request.urlopen(url) as response, open(output_path, 'wb') as out_file:
        shutil.copyfileobj(response, out_file)


def download_and_extract_gz(url, output_path):
    """Download a gzipped file and extract it."""
    logger.info(f"Downloading and extracting {url} to {output_path}")
    temp_file = tempfile.NamedTemporaryFile(delete=False).name
    try:
        download_file(url, temp_file)
        with gzip.open(temp_file, 'rb') as f_in, open(output_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    finally:
        if os.path.exists(temp_file):
            os.unlink(temp_file)


def scenario_1_genome_annotation():
    """
    Scenario 1: Genome Annotation Pipeline
    
    This scenario demonstrates how to use diamondonpy for annotating coding sequences (CDS)
    from a newly sequenced genome against a reference protein database (like UniProt).
    
    Workflow:
    1. Download a small reference protein database (subset of UniProt)
    2. Create a DIAMOND database from the reference proteins
    3. Simulate a newly sequenced genome with coding sequences
    4. Use blastx to translate the nucleotide sequences to proteins and align against reference
    5. Filter the results to find the best hits
    6. Generate a simple annotation report
    """
    logger.info("Running Scenario 1: Genome Annotation Pipeline")
    
    # Initialize Diamond
    diamond = Diamond()
    
    # Step 1: Download reference protein database (using a small example)
    uniprot_sample = DATA_DIR / "uniprot_sample.fasta"
    if not uniprot_sample.exists():
        # Note: In a real scenario, you would use the full UniProt database
        # This is a small sample for demonstration purposes
        # download_file(
        #     "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%28proteome%3AUP000464024%29%29",
        #     uniprot_sample
        # )
        # In a real scenario, you would add some actual protein sequences here
        with open(uniprot_sample, 'a') as f:
            f.write(">sp|P0DTC2|SPIKE_SARS2 Spike glycoprotein OS=Severe acute respiratory syndrome coronavirus 2 OX=2697049 GN=S PE=1 SV=1\n")
            f.write("MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHV\n")
            f.write("SGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCND\n")
    
    # Step 2: Create a DIAMOND database
    db_path = DATA_DIR / "uniprot_sample.dmnd"
    if not db_path.exists():
        diamond.makedb(
            db=str(db_path),
            input_file=str(uniprot_sample),
            threads=4
        )
        logger.info(f"Created DIAMOND database at {db_path}")
    
    # Step 3: Create a sample of nucleotide sequences (simulating newly sequenced genome)
    new_genome_cds = DATA_DIR / "new_genome_cds.fasta"
    if not new_genome_cds.exists():
        with open(new_genome_cds, 'w') as f:
            # Example CDS sequence encoding a protein similar to spike
            f.write(">contig1_cds1 putative viral protein\n")
            f.write("ATGTTCGTGTTCCTGGTGCTGCTGCCACTGGTGTCCAGCCAGTGCGTGAACCTGACCACGCGGACCCAGCTGCCACCGGCGTACACCAACTCCTTCACGCGGGGCGTCTACTACCCTGACAAGGTGTTCCGGAGCTCGGTGCTGCACAGCACCCAGGACCTGTTCCTGCCCTTCTTCAGCAACGTGACCTGGTTCCACGCCATCCACGTGAGCGGCACCAACGGCACCAAGCGGTTCGACAACCCCGTGCTGCCCTTCAACGACGGCGTGTACTTCGCCTCGACGGAGAAGAGCAACATCATCCGCGGCTGGATCTTCGGCACCACACTGGACTCCAAGACCCAGTCCCTGCTGATCGTGAACAACGCCACCAACGTGGTGATCAAGGTGTGCGAGTTCCAGTTCTGCAACGAC\n")
    
    # Step 4: Run BLASTX to translate and align the sequences
    logger.info("Running BLASTX to align translated CDS against protein database")
    results = diamond.blastx(
        db=str(db_path),
        query=str(new_genome_cds),
        evalue=1e-5,
        max_target_seqs=5,
        outfmt="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe stitle",
        threads=4
    )
    
    # Step 5: Filter and analyze results
    logger.info(f"Found {len(results)} hits")
    if len(results) > 0:
        # Sort by bitscore (descending)
        results = results.sort_values(by='bitscore', ascending=False)
        
        # Display top hits
        logger.info("Top hits for annotation:")
        print(results[['qseqid', 'sseqid', 'pident', 'evalue', 'bitscore', 'stitle']].head())
        
        # Generate an annotation report
        report_path = DATA_DIR / "annotation_report.csv"
        results.to_csv(report_path, index=False)
        logger.info(f"Annotation report saved to {report_path}")
    else:
        logger.warning("No hits found for the sample sequences")
    
    logger.info("Genome annotation pipeline completed")


def scenario_2_protein_clustering():
    """
    Scenario 2: Protein Family Analysis through Clustering
    
    This scenario demonstrates how to use diamondonpy for clustering similar proteins
    into protein families or groups based on sequence similarity.
    
    Workflow:
    1. Create/download a dataset of related proteins (e.g., kinases, transporters)
    2. Use clustering to group similar proteins
    3. Analyze the clusters to identify protein families
    4. Visualize cluster sizes and relationships
    """
    logger.info("Running Scenario 2: Protein Family Analysis through Clustering")
    
    # Initialize Diamond
    diamond = Diamond()
    
    # Step 1: Create a sample protein dataset
    protein_dataset = DATA_DIR / "protein_families.fasta"
    if not protein_dataset.exists():
        with open(protein_dataset, 'w') as f:
            # Add some sample protein sequences (in a real-world scenario, you'd have many more)
            # Family 1: Transport proteins
            f.write(">tr|A0A1X7SBP1|A0A1X7SBP1_9SPHI Transporter OS=Sphingobacterium sp.\n")
            f.write("MKKLAILAATLVGLGLAACSQAAQAETYTVKLGNDHIKADFRVGQPAELVEGKDLSGLGLQAGDHIKAGF\n")
            f.write("EVGSAENVAAARDALASQLGVSYAQVKEAITSAPAAVAAHGALAKTEMKDTK\n")
            
            f.write(">tr|A0A1C4KZ81|A0A1C4KZ81_9ACTN Transporter OS=Streptomyces sp.\n")
            f.write("MKTQRNLLSGLLLVAVLVLGGCGSEQAATANASLSLEGRDYKADARVKTGPDTISEGKEIDGLGLTPGAYAKAGF\n")
            f.write("DVKSALQLSEAIKALAAAPSANYKDVKALTESAPAAMKAYGILAKTAIK\n")
            
            # Family 2: Kinases
            f.write(">sp|P00531|SRC_CHICK Proto-oncogene tyrosine-protein kinase Src OS=Gallus gallus\n")
            f.write("MGSNKSKPKDASQRRRSLEPAENVHGAGGGAFPASQTPSKPASADGHRGPSAAFAPAAAEPKLFGGFNSSDTVTS\n")
            f.write("PGQAEGLCPPVPKLAKDKGSAQEPQPTPTYQALLQPHIVYCSKTGDYWEIPRESLRLEVKLGQGCFGEVWMGTW\n")
            
            f.write(">sp|P12931|SRC_HUMAN Proto-oncogene tyrosine-protein kinase Src OS=Homo sapiens\n")
            f.write("MGSNKSKPKDASQRRRSLEPAENVHGAGGGAFPASQTPSKPASADGHRGPSAAFAPAAAEPKLFGGFNSSDTVTS\n")
            f.write("PGQAEPVQPTPTYQALLQPHIVYCSKTDGYWEIPRESLRLEVKLGQGCFGEVWMGTWNGTTRVAIKTLKPGTMSP\n")
            
            # Family 3: Enzymes
            f.write(">sp|P00762|TRY1_RAT Anionic trypsin-1 OS=Rattus norvegicus\n")
            f.write("MKTFIFLALLGAAVAFPLEDDDKIVGGYTCGANTVPYQVSLNSGYHFCGGSLINSQWVVSAAHCYKSGIQVRLGE\n")
            f.write("DNINVVEGNEQFISASKSIVHPSYNSNTLNNDIMLIKLKSAASLNSRVASISLPTSCA\n")
            
            f.write(">sp|P00766|CTRA_BOVIN Chymotrypsinogen A OS=Bos taurus\n")
            f.write("CGVPAIQPVLSGLSRIVNGEEAVPGSWPWQVSLQDKTGFHFCGGSLINENWVVTAAHCGVTTSDVVVAGEFDQGS\n")
            f.write("DKEKIQKLKIAKVFKNSKYNSLTINNDITLLKLSTAASFSQTVSAVCLPSASDDFAAGTTCVTTGWGLTRYTNANT\n")
    
    # Step 2: Run clustering
    logger.info("Clustering proteins by sequence similarity")
    db_path = DATA_DIR / "protein_families.dmnd"
    
    # Create database
    diamond.makedb(
        db=str(db_path),
        input_file=str(protein_dataset),
        threads=4
    )
    
    # Perform clustering
    clusters = diamond.cluster(
        db=str(db_path),
        threads=4,
    )
    
    # Step 3: Analyze the clusters
    logger.info(f"Found {len(clusters)} protein entries in {clusters['cluster_number'].nunique()} clusters")
    
    # Count proteins per cluster
    cluster_sizes = clusters.groupby('cluster_number').size().reset_index(name='size')
    print(clusters)
    # Print cluster information
    print("\nCluster information:")
    for cluster_num, group in clusters.groupby('cluster_number'):
        print(f"Cluster {cluster_num}: {len(group)} proteins")
        print(f"  Representative: {group.iloc[0]['representative']}")
        members = group['sequence_id'].tolist()
        print(f"  Members: {', '.join(members[:3])}" + ("..." if len(members) > 3 else ""))
    
    # Step 4: Visualize cluster sizes
    if 'matplotlib.pyplot' in sys.modules:
        plt.figure(figsize=(10, 6))
        plt.bar(cluster_sizes['cluster_number'].astype(str), cluster_sizes['size'])
        plt.xlabel('Cluster')
        plt.ylabel('Number of proteins')
        plt.title('Protein Family Cluster Sizes')
        plt.tight_layout()
        
        # Save figure
        plt.savefig(DATA_DIR / "cluster_sizes.png")
        logger.info(f"Cluster size visualization saved to {DATA_DIR / 'cluster_sizes.png'}")
    
    logger.info("Protein clustering analysis completed")


def scenario_3_ortholog_identification():
    """
    Scenario 3: Ortholog Identification Across Species
    
    This scenario demonstrates how to use the bidirectional best hit approach
    to identify orthologous genes between two species.
    
    Workflow:
    1. Create sample protein datasets for two species
    2. Create DIAMOND databases for both species
    3. Perform bidirectional best hit analysis
    4. Analyze and visualize the results
    """
    logger.info("Running Scenario 3: Ortholog Identification")
    
    # Initialize Diamond
    diamond = Diamond()
    
    # Step 1: Create sample datasets for two species
    species1_proteins = DATA_DIR / "species1_proteins.fasta"
    species2_proteins = DATA_DIR / "species2_proteins.fasta"
    
    if not species1_proteins.exists():
        with open(species1_proteins, 'w') as f:
            # Human proteins
            f.write(">sp|P68871|HBB_HUMAN Hemoglobin subunit beta OS=Homo sapiens\n")
            f.write("MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSD\n")
            f.write("GLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH\n")
            
            f.write(">sp|P01308|INS_HUMAN Insulin OS=Homo sapiens\n")
            f.write("MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPG\n")
            f.write("AGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN\n")
            
            f.write(">sp|P35968|VGFR2_HUMAN Vascular endothelial growth factor receptor 2 OS=Homo sapiens\n")
            f.write("MQSKVLLAVALWLCVETRAASVGLPSVSLDLPRLSIQKDILTIKANTTLQITCRGQRDLDWLWPNNQSGSEQR\n")
            f.write("VEVTECSDGLFCKTLTIPKVIGNDTGAYKCFYRETDLASVIYVYVQDYRSPFIASVSDQHGVVYITENKNKTVV\n")
    
    if not species2_proteins.exists():
        with open(species2_proteins, 'w') as f:
            # Mouse proteins (orthologs to the human proteins)
            f.write(">sp|P02088|HBB1_MOUSE Hemoglobin subunit beta-1 OS=Mus musculus\n")
            f.write("MVHLTDAEKAAVSCLWGKVNSDEVGGEALGRLLVVYPWTQRYFDSFGDLSSASAIMGNAKVKAHGKKVITAFND\n")
            f.write("GLNHLDSLKGTFASLSELHCDKLHVDPENFRLLGNMIVIVLGHHLGKDFTPAAQAAFQKVVAGVATALAHKYH\n")
            
            f.write(">sp|P01326|INS2_MOUSE Insulin-2 OS=Mus musculus\n")
            f.write("MALWMRFLPLLALLVLWEPKPAQAFVKQHLCGPHLVEALYLVCGERGFFYTPMSRREVEDPQVAQLELGGGPG\n")
            f.write("AGDLQTLALEVAQQKRGIVDQCCTSICSLYQLENYCN\n")
            
            f.write(">sp|P35918|VGFR2_MOUSE Vascular endothelial growth factor receptor 2 OS=Mus musculus\n")
            f.write("MQSKVLLAAALWFCMETQAASVGLPAVSLDLPRLSIQRSVLITRGNSTTLEIQCPGQRDLDWLWPNNQSGSER\n")
            f.write("VEVLECLDGLICTRLNIPRITGNDTGTYRCTFRETDIMASIYVYVKNYGSPLIASISDKHGVVYITENRKATIV\n")
            
            # An additional mouse protein without a direct ortholog in our human set
            f.write(">sp|P11404|FABPH_MOUSE Fatty acid-binding protein, heart OS=Mus musculus\n")
            f.write("MVDAFLGTWKLVDSKNFDDYMKSLGVGFATRQVASMRTVVPTQEKVQQKMTITFKNVVNFEFCAGKLSVHIQKG\n")
            f.write("DDHDMFKKMVGLAKRKLPCSVTKSIQVGQITVADLKGKTTVTVRELRGQVLVNVQFGEENTQMTFKGIKSVTEL\n")
    
    # Step 3: Perform bidirectional best hit analysis
    logger.info("Performing bidirectional best hit analysis to identify orthologs")
    
    orthologs = diamond.bidirectional_best_hit(
        fasta1=str(species1_proteins),
        fasta2=str(species2_proteins),
        evalue=1e-5,
        threads=4
    )
    
    # Step 4: Analyze results
    logger.info(f"Found {len(orthologs)} potential ortholog pairs")
    
    # Display ortholog pairs
    print("\nOrtholog pairs identified:")
    for idx, row in orthologs.iterrows():
        print(f"Species 1: {row['qseqid']} ↔ Species 2: {row['sseqid']} (E-value: {row['evalue']:.2e}, Identity: {row['pident']:.1f}%)")
    
    # Save results
    orthologs.to_csv(DATA_DIR / "ortholog_pairs.csv", index=False)
    logger.info(f"Ortholog pair data saved to {DATA_DIR / 'ortholog_pairs.csv'}")
    
    logger.info("Ortholog identification completed")


def scenario_4_metagenome_analysis():
    """
    Scenario 4: Metagenomic Analysis
    
    This scenario demonstrates how to analyze metagenomic data to identify 
    the taxonomic and functional composition of a microbial community.
    
    Workflow:
    1. Create/download a sample metagenomic dataset (contigs or reads)
    2. Use a reference protein database 
    3. Translate and align the metagenomic sequences using blastx
    4. Analyze the taxonomic and functional assignments
    5. Generate summary statistics and visualizations
    """
    logger.info("Running Scenario 4: Metagenomic Analysis")
    
    # Initialize Diamond
    diamond = Diamond()
    
    # Step 1: Create a sample metagenomic dataset
    metagenome_file = DATA_DIR / "metagenome_sample.fasta"
    if not metagenome_file.exists():
        with open(metagenome_file, 'w') as f:
            # Simulated metagenomic contigs (in a real-world scenario, you'd have thousands)
            f.write(">contig1 length=1200\n")
            f.write("ATGAAGAAACTTATCGCACTCGCAGCAACACTCGTTGGCCTCGGCCTTGCAGCCTGCTCACAACCTTCACAAGCC\n")
            f.write("GCTCAAGCGGCGGAAACCTATACCGTCAAACTGGGCAACGACCACATCAAAGCCGACTTCCGCGTCGGCCAGCCT\n")
            f.write("GCGGAACTCGTCGAAGGCAAGGACCTGTCCGGCCTCGGCCTGCAGGCAGGCGACCACATCAAGGCAGGCTTCGAA\n")
            f.write("GTCGGCTCCGCCGAAAACGTCGCAGCCGCACGCGACGCACTCGCATCGCAGCTCGGCGTCAGCTACGCACAGGTC\n")
            
            f.write(">contig2 length=900\n")
            f.write("ATGAAAACACAGCGCAATCTGCTGTCCGGCCTGCTGCTCGTCGCAGTGCTCGTCCTCGGCGGCTGCGGCAGCGAA\n")
            f.write("CAAGCGGCGACCGCAAACGCAAGCCTGAGCCTCGAGGGCCGCGACTACAAGGCCGACGCACGCGTCAAGACCGGC\n")
            f.write("CCCGACACCATCTCCGAGGGCAAGGAGATCGACGGCCTCGGCCTGACACCCGGCGCATACGCCAAGGCAGGCTTC\n")
            
            f.write(">contig3 length=1500\n")
            f.write("ATGGGCTCAAACAAGTCCAAGCCCAAGGACGCCTCCCAGCGCCGCCGCTCCCTCGAGCCCGCCGAGAACGTCCAC\n")
            f.write("GGCGCCGGCGGCGGCGCCTTCCCCGCCTCCCAGACGCCCTCCAAGCCCGCCTCCGCCGACGGCCACCGCGGCCCG\n")
            f.write("TCCGCGGCCTTCGCCCCCGCCGCCGCCGAGCCCAAGCTCTTCGGCGGCTTCAACTCCTCCGACACCGTCACCTCC\n")
    
    # Step 2: Create or use a reference database
    # For simplicity, we'll reuse the protein database from scenario 1
    uniprot_sample = DATA_DIR / "uniprot_sample.fasta"
    if not uniprot_sample.exists():
        logger.info("Please run scenario 1 first to create the reference database")
        return
    
    db_path = DATA_DIR / "uniprot_sample.dmnd"
    if not db_path.exists():
        logger.info("Please run scenario 1 first to create the DIAMOND database")
        return
    
    # Step 3: Translate and align metagenomic sequences
    logger.info("Translating and aligning metagenomic sequences against reference database")
    results = diamond.blastx(
        db=str(db_path),
        query=str(metagenome_file),
        evalue=1e-3,  # Less stringent for metagenomic data
        max_target_seqs=5,
        outfmt="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe stitle",
        threads=4
    )
    
    # Step 4: Analyze the results
    if len(results) > 0:
        logger.info(f"Found {len(results)} hits for {results['qseqid'].nunique()} contigs")
        
        # Group by contig and get best hit
        best_hits = results.loc[results.groupby('qseqid')['bitscore'].idxmax()]
        
        # Extract taxonomic information (in a real scenario, you would parse this from the full UniProt headers)
        # Here we'll just simulate it with the limited data we have
        best_hits['organism'] = best_hits['stitle'].str.extract(r'OS=(.*?)(?:\s+OX=|\s+GN=|\s+PE=|\s+SV=|$)')
        
        # Count hits per organism
        org_counts = best_hits['organism'].value_counts().reset_index()
        org_counts.columns = ['Organism', 'Count']
        
        # Display taxonomic summary
        print("\nTaxonomic composition:")
        print(org_counts)
        
        # Save results
        results.to_csv(DATA_DIR / "metagenome_analysis.csv", index=False)
        logger.info(f"Metagenomic analysis results saved to {DATA_DIR / 'metagenome_analysis.csv'}")
        
        # Visualize (in a real scenario, you'd create more detailed visualizations)
        if 'matplotlib.pyplot' in sys.modules and len(org_counts) > 0:
            plt.figure(figsize=(10, 6))
            plt.bar(org_counts['Organism'], org_counts['Count'])
            plt.xlabel('Organism')
            plt.ylabel('Number of contigs')
            plt.title('Taxonomic Composition of Metagenome')
            plt.xticks(rotation=45, ha='right')
            plt.tight_layout()
            
            # Save figure
            plt.savefig(DATA_DIR / "taxonomic_composition.png")
            logger.info(f"Taxonomic visualization saved to {DATA_DIR / 'taxonomic_composition.png'}")
    else:
        logger.warning("No hits found for the metagenomic sequences")
    
    logger.info("Metagenomic analysis completed")


def scenario_5_domain_analysis():
    """
    Scenario 5: Protein Domain Analysis
    
    This scenario demonstrates how to identify conserved domains in protein sequences
    by aligning them against a database of known domains and analyzing the alignment results.
    
    Workflow:
    1. Create a set of protein sequences
    2. Use a domain database (or simulate one)
    3. Use blastp to align proteins against domain database
    4. Analyze and visualize domain architecture
    """
    logger.info("Running Scenario 5: Protein Domain Analysis")
    
    # Initialize Diamond
    diamond = Diamond()
    
    # Step 1: Create a sample protein dataset
    protein_file = DATA_DIR / "query_proteins.fasta"
    if not protein_file.exists():
        with open(protein_file, 'w') as f:
            # Sample protein with multiple domains
            f.write(">protein1 multi-domain protein\n")
            f.write("MGSNKSKPKDASQRRRSLEPAENVHGAGGGAFPASQTPSKPASADGHRGPSAAFAPAAAEPKLFGGFNSSDTVTS\n")
            f.write("PGQAEPVQPTPTYQALLQPHIVYCSKTDGYWEIPRESLRLEVKLGQGCFGEVWMGTWNGTTRVAIKTLKPGTMSP\n")
            f.write("EAFLQEAQIMKKLRHEKLVQLYAVVSEEPIYIVTEYMSKGSLLDFLKGETGKYLRLPQLVDMAAQIASGMAYVER\n")
            f.write("MNYALQTGGGRALSEEDKFYWASEPGNPQKLSQSQPVGQPVPRVYSTTAVQPQHLSGYRPERPGAADETHGNVT\n")
            
            # Another protein with different domains
            f.write(">protein2 signaling protein\n")
            f.write("MAAAAAAGAGPEMVRGQVFDVGPRYTNLSYIGEGAYGMVCSAYDNVNKVRVAIKKISPFEHQTYCQRTLREIKIL\n")
            f.write("LRFRHENIIGINDIIRAPTIEQMKDVYIVQDLMETDLYKLLKTQHLSNDHICYFLYQILRGLKYIHSANVLHRDL\n")
            f.write("KPSNLLLNTTCDLKICDFGLARVADPDHDHTGFLTEYVATRWYRAPEIMLNSKGYTKSIDIWSVGCILAEMLSNR\n")
            
            # A simpler protein
            f.write(">protein3 small protein\n")
            f.write("MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVA\n")
            f.write("PAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDS\n")
    
    # Step 2: Create a simple domain database
    domain_file = DATA_DIR / "domain_db.fasta"
    if not domain_file.exists():
        with open(domain_file, 'w') as f:
            # Some example domains
            f.write(">domain1 SH2 domain\n")
            f.write("PQTPTYQALLQPHIVYCSKTDGYWEIPRESLRLEVKLGQGCFG\n")
            
            f.write(">domain2 Protein kinase domain\n")
            f.write("EVWMGTWNGTTRVAIKTLKPGTMSPEAFLQEAQIMKKLRHEKLVQLYAVVS\n")
            
            f.write(">domain3 MAPK domain\n")
            f.write("GAYGMVCSAYDNVNKVRVAIKKISPFEHQTYCQRTLREIKILLRFRHEN\n")
            
            f.write(">domain4 p53 domain\n")
            f.write("PAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTY\n")
    
    # Create database
    domain_db = DATA_DIR / "domain_db.dmnd"
    if not domain_db.exists():
        diamond.makedb(
            db=str(domain_db),
            input_file=str(domain_file),
            threads=4
        )
    
    # Step 3: Run blastp to find domains
    logger.info("Aligning proteins against domain database")
    results = diamond.blastp(
        db=str(domain_db),
        query=str(protein_file),
        evalue=1e-2,
        outfmt="6 qseqid sseqid pident length qstart qend sstart send evalue bitscore stitle",
        threads=4
    )
    
    # Step 4: Analyze domain hits
    if len(results) > 0:
        logger.info(f"Found {len(results)} domain hits for {results['qseqid'].nunique()} proteins")
        
        # Sort by protein and then by start position
        results = results.sort_values(by=['qseqid', 'qstart'])
        
        # Display domain architecture for each protein
        print("\nDomain architecture:")
        for protein, group in results.groupby('qseqid'):
            print(f"\n{protein}:")
            for _, hit in group.iterrows():
                domain = hit['stitle'].split()[0] if pd.notna(hit['stitle']) else hit['sseqid']
                print(f"  {domain} at positions {hit['qstart']}-{hit['qend']} ({hit['length']} aa, {hit['pident']:.1f}% identity)")
        
        # Save results
        results.to_csv(DATA_DIR / "domain_analysis.csv", index=False)
        logger.info(f"Domain analysis results saved to {DATA_DIR / 'domain_analysis.csv'}")
        
        # Visualize domain architecture (simple text-based visualization)
        print("\nDomain visualization:")
        for protein, group in results.groupby('qseqid'):
            # Get protein length (in a real scenario, you'd get this from FASTA)
            protein_length = group['qend'].max()
            
            # Create a simple text-based visualization
            visualization = ['-'] * protein_length
            for _, hit in group.iterrows():
                domain_id = ''.join(c[0] for c in hit['sseqid'].split('_'))
                start, end = int(hit['qstart']), int(hit['qend'])
                for i in range(start, min(end + 1, protein_length + 1)):
                    visualization[i-1] = domain_id
            
            # Print the visualization
            print(f"\n{protein}:")
            print(''.join(visualization))
    else:
        logger.warning("No domain hits found")
    
    logger.info("Protein domain analysis completed")


def main():
    """Main function to run the scenarios."""
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Run real-world scenarios for diamondonpy")
    parser.add_argument('--scenario', type=int, choices=range(1, 6), default=1,
                        help="Scenario to run (1-5)")
    parser.add_argument('--all', action='store_true', help="Run all scenarios")
    args = parser.parse_args()
    
    # Dictionary mapping scenario numbers to functions
    scenarios = {
        1: scenario_1_genome_annotation,
        2: scenario_2_protein_clustering,
        3: scenario_3_ortholog_identification,
        4: scenario_4_metagenome_analysis,
        5: scenario_5_domain_analysis
    }
    
    if args.all:
        # Run all scenarios
        for scenario_func in scenarios.values():
            scenario_func()
            print("\n" + "="*80 + "\n")
    else:
        # Run the specified scenario
        scenarios[args.scenario]()


if __name__ == "__main__":
    main() 