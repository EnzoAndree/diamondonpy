"""Integration tests for the Diamond package."""
import pytest
import os
import pandas as pd
from diamondonpy import Diamond, DiamondError

@pytest.fixture
def diamond():
    """Return a Diamond instance for testing."""
    try:
        d = Diamond()
        return d
    except DiamondError:
        pytest.skip("DIAMOND not installed or not accessible")

@pytest.fixture
def test_data(tmp_path):
    """Create test FASTA files and paths."""
    # Create a query FASTA file with DNA sequences for blastx
    query_fasta = tmp_path / "query.fasta"
    query_fasta.write_text(
        ">seq1\nATGAAACTGGTGGTGGTGCCGGAGACGCGGCCGAACCCGAACGGCTACAAGTTCTCGTTCAAGAAGGTGAAGGAGGTGCTGAAGTCGCTGCCGGAGGAGAAGCGCAAGGCGTACGAGGAGCTGGCGCGGGAGCTGGGCCTGAACCCGGAGGAGGTGGCGCGGCGGCTGAAGGCGAAGCTGGAGGAGCTGGGCCTG\n"
        ">seq2\nATGATCGCGGTGGCGGTGGCGCTGGCGGGCGGCGCGACGCAGGCGTTCGCGAACCCGCTGCCGGGCCAGTTCGGCAAGACGCTGACGCTGAAGGGCAAGGCGAAGCTGGTGATCGGCGACGAGGTGCGGCTGACGAACCCGCTGGGCGTGCCGTCGCGGATGGAGCGGGTGCGGCTGGTGGACCTGTCGAACCCGTAC\n"
    )
    
    # Create a reference FASTA file with protein sequences
    ref_fasta = tmp_path / "ref.fasta"
    ref_fasta.write_text(
        ">ref1\nMKLVVVPETRPNPNGYKFSFKKVKEVLKSLPEEKRKAYEELARELGLNPEEVARRLKAKLEELGL\n"
        ">ref2\nMIAVAVALAGGATQAFANPLPGQFGKTLTLKGKAKLVIGDEVRLTNPLGVPSRMERVRLVDLSNPY\n"
        ">ref3\nMKLIIIPDTRPNPNGYKFSFKKVKEVLKSLPEEKRKAYEELARELGLNPEEVARRLKAKLEELGL\n"
    )
    
    # Create paths for output files
    db_path = tmp_path / "ref.dmnd"
    out_path = tmp_path / "results.txt"
    daa_path = tmp_path / "alignment.daa"
    
    # Add FASTA content as strings
    fasta_strings = {
        'query_dna': ">seq1\nATGAAACTGGTGGTGGTGCCGGAGACGCGGCCGAACCCGAACGGCTACAAGTTCTCGTTCAAGAAGGTGAAGGAGGTGCTGAAGTCGCTGCCGGAGGAGAAGCGCAAGGCGTACGAGGAGCTGGCGCGGGAGCTGGGCCTGAACCCGGAGGAGGTGGCGCGGCGGCTGAAGGCGAAGCTGGAGGAGCTGGGCCTG\n",
        'ref_protein': ">ref1\nMKLVVVPETRPNPNGYKFSFKKVKEVLKSLPEEKRKAYEELARELGLNPEEVARRLKAKLEELGL\n"
    }
    
    return {
        'query': str(query_fasta),
        'ref': str(ref_fasta),
        'db': str(db_path),
        'out': str(out_path),
        'daa': str(daa_path),
        'fasta_strings': fasta_strings
    }

@pytest.mark.integration
def test_makedb(diamond, test_data):
    """Test database creation with real FASTA file."""
    output = diamond.makedb(
        db=test_data['db'],
        input_file=test_data['ref']
    )
    assert os.path.exists(test_data['db'])

@pytest.mark.integration
def test_blastp(diamond, test_data):
    """Test protein BLAST search with real sequences."""
    # First create the database
    diamond.makedb(db=test_data['db'], input_file=test_data['ref'])
    
    # Run BLASTP search with protein sequences
    df = diamond.blastp(
        db=test_data['db'],
        query=test_data['ref'],  # Use protein sequences for blastp
        evalue=1e-10
    )
    
    assert isinstance(df, pd.DataFrame)
    assert len(df) > 0  # Should find at least one hit
    assert 'qseqid' in df.columns
    assert 'sseqid' in df.columns
    assert 'pident' in df.columns
    
    # Check for exact match
    exact_match = df[df['pident'] == 100.0]
    assert not exact_match.empty, "Should find at least one 100% identity match"

@pytest.mark.integration
def test_blastx(diamond, test_data):
    """Test translated BLAST search with real sequences."""
    # First create the database
    diamond.makedb(db=test_data['db'], input_file=test_data['ref'])
    
    # Run BLASTX search
    df = diamond.blastx(
        db=test_data['db'],
        query=test_data['query'],
        evalue=1e-10
    )
    
    assert isinstance(df, pd.DataFrame)
    assert len(df) > 0
    assert 'qseqid' in df.columns
    assert 'sseqid' in df.columns
    assert 'pident' in df.columns

@pytest.mark.integration
def test_output_formats(diamond, test_data):
    """Test different output formats."""
    # Create database
    diamond.makedb(db=test_data['db'], input_file=test_data['ref'])
    
    # Test tabular output with custom format
    df = diamond.blastp(
        db=test_data['db'],
        query=test_data['ref'],  # Use protein query for blastp
        outfmt=6,
        evalue=1e-10
    )
    
    assert isinstance(df, pd.DataFrame)
    assert len(df) > 0
    assert 'qseqid' in df.columns
    
    # Test output to file
    diamond.blastp(
        db=test_data['db'],
        query=test_data['ref'],
        out=test_data['out'],
        evalue=1e-10
    )
    assert os.path.exists(test_data['out'])

@pytest.mark.integration
def test_clustering(diamond, test_data):
    """Test sequence clustering."""
    # Create database
    diamond.makedb(db=test_data['db'], input_file=test_data['ref'])
    
    # Test clustering
    df = diamond.cluster(
        db=test_data['db'],
        approx_id=90.0
    )
    
    assert isinstance(df, pd.DataFrame)
    assert len(df) > 0
    assert 'representative' in df.columns
    assert 'sequence_id' in df.columns
    assert 'cluster_number' in df.columns
    assert 'cluster_size' in df.columns
    
    # Test linclust
    df = diamond.linclust(
        db=test_data['db'],
        approx_id=90.0
    )
    
    assert isinstance(df, pd.DataFrame)
    assert len(df) > 0
    assert 'representative' in df.columns
    assert 'sequence_id' in df.columns
    assert 'cluster_number' in df.columns
    assert 'cluster_size' in df.columns

@pytest.mark.integration
def test_makedb_from_string(diamond, test_data):
    """Test database creation with FASTA string input."""
    output = diamond.makedb(
        db=test_data['db'],
        input_file=test_data['fasta_strings']['ref_protein']
    )
    assert os.path.exists(test_data['db'])

@pytest.mark.integration
def test_blastp_with_string_input(diamond, test_data):
    """Test protein BLAST search with FASTA string input."""
    # First create the database
    diamond.makedb(db=test_data['db'], input_file=test_data['ref'])
    
    # Run BLASTP search with protein sequence string
    df = diamond.blastp(
        db=test_data['db'],
        query=test_data['fasta_strings']['ref_protein'],
        evalue=1e-10
    )
    
    assert isinstance(df, pd.DataFrame)
    assert len(df) > 0
    assert 'qseqid' in df.columns
    assert 'sseqid' in df.columns
    assert 'pident' in df.columns
    
    # Check for exact match
    exact_match = df[df['pident'] == 100.0]
    assert not exact_match.empty, "Should find at least one 100% identity match"

@pytest.mark.integration
def test_blastx_with_string_input(diamond, test_data):
    """Test translated BLAST search with FASTA string input."""
    # First create the database
    diamond.makedb(db=test_data['db'], input_file=test_data['ref'])
    
    # Run BLASTX search with DNA sequence string
    df = diamond.blastx(
        db=test_data['db'],
        query=test_data['fasta_strings']['query_dna'],
        evalue=1e-10
    )
    
    assert isinstance(df, pd.DataFrame)
    assert len(df) > 0
    assert 'qseqid' in df.columns
    assert 'sseqid' in df.columns
    assert 'pident' in df.columns

@pytest.mark.integration
def test_mixed_input_types(diamond, test_data):
    """Test using both file paths and string inputs in the same workflow."""
    # Create database from string input
    diamond.makedb(
        db=test_data['db'],
        input_file=test_data['fasta_strings']['ref_protein']
    )
    
    # Run BLASTP with file input against string-created database
    df_file = diamond.blastp(
        db=test_data['db'],
        query=test_data['ref'],
        evalue=1e-10
    )
    
    # Run BLASTP with string input against same database
    df_string = diamond.blastp(
        db=test_data['db'],
        query=test_data['fasta_strings']['ref_protein'],
        evalue=1e-10
    )
    
    assert isinstance(df_file, pd.DataFrame)
    assert isinstance(df_string, pd.DataFrame)
    assert len(df_file) > 0
    assert len(df_string) > 0
    
    # Both should find exact matches
    assert not df_file[df_file['pident'] == 100.0].empty
    assert not df_string[df_string['pident'] == 100.0].empty 