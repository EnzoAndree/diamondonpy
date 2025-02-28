"""Unit tests for the Diamond class."""
import pytest
import pandas as pd
import os
from diamondonpy import Diamond, DiamondError
import subprocess

@pytest.mark.unit
def test_diamond_init(mock_successful_run):
    """Test Diamond class initialization."""
    diamond = Diamond()
    assert diamond.executable == "diamond"
    
    # Test custom executable path
    diamond = Diamond(executable="/custom/path/diamond")
    assert diamond.executable == "/custom/path/diamond"

@pytest.mark.unit
def test_verify_installation(mocker):
    """Test installation verification."""
    # Should not raise with successful run
    mock_run = mocker.patch('subprocess.run')
    mock_run.return_value.stdout = "diamond version 2.1.0"
    mock_run.return_value.stderr = ""
    mock_run.return_value.returncode = 0
    Diamond()
    
    # Should raise with failed run
    mock_run.side_effect = subprocess.CalledProcessError(1, ["diamond", "--version"], stderr="mock error")
    with pytest.raises(DiamondError):
        Diamond()

@pytest.mark.unit
def test_build_command(diamond):
    """Test command building."""
    cmd = diamond._build_command("makedb", {
        "db": "test.dmnd",
        "input_file": "test.fasta",
        "threads": 4
    })
    
    assert cmd == [
        "diamond",
        "makedb",
        "--db", "test.dmnd",
        "--in", "test.fasta",
        "--threads", "4"
    ]
    
    # Test boolean flag
    cmd = diamond._build_command("view", {
        "daa": "test.daa",
        "verbose": True
    })
    
    assert cmd == [
        "diamond",
        "view",
        "--daa", "test.daa",
        "--verbose"
    ]
    
    # Test outfmt validation
    with pytest.raises(ValueError, match="Invalid column names"):
        diamond._build_command("blastp", {
            "db": "test.dmnd",
            "query": "test.fasta",
            "outfmt": "6 invalid_column"
        })

@pytest.mark.unit
def test_parse_blast_output(diamond):
    """Test parsing of BLAST tabular output."""
    # Test default format
    mock_output = "seq1\tref1\t98.5\t100\t2\t0\t1\t100\t1\t100\t1e-50\t200.0\n"
    df = diamond._parse_blast_output(mock_output)
    
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 1
    assert list(df.columns) == [
        'qseqid', 'sseqid', 'pident', 'length', 'mismatch',
        'gapopen', 'qstart', 'qend', 'sstart', 'send',
        'evalue', 'bitscore'
    ]
    assert df.iloc[0]['qseqid'] == 'seq1'
    assert df.iloc[0]['pident'] == 98.5
    
    # Test custom column format
    mock_output = "seq1\tref1\t98.5\n"
    df = diamond._parse_blast_output(mock_output, outfmt="6 qseqid sseqid pident")
    
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 1
    assert list(df.columns) == ['qseqid', 'sseqid', 'pident']
    assert df.iloc[0]['qseqid'] == 'seq1'
    assert df.iloc[0]['pident'] == 98.5
    
    # Test empty output
    mock_output = ""
    df = diamond._parse_blast_output(mock_output)
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 0
    
    # Test invalid format - providing just "6" (as a string)
    # This is now handled as the default format of 6 columns
    with pytest.raises(ValueError) as excinfo:
        diamond._parse_blast_output("", outfmt="invalid")
    assert "Invalid output format specification" in str(excinfo.value)
    
    # Test parsing error for mismatched columns
    # This should now be handled by the flexible parser
    mock_output = "invalid\tformat\twith\textra\tcolumns\n"
    df = diamond._parse_blast_output(mock_output)
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 1

@pytest.mark.unit
def test_makedb(diamond, mock_successful_run, test_files):
    """Test makedb command."""
    output = diamond.makedb(
        db=test_files['db'],
        input_file=test_files['fasta'],
        threads=4
    )
    
    assert isinstance(output, str)
    assert output == "mock output"
    assert mock_successful_run.call_count >= 1
    cmd = mock_successful_run.call_args_list[-1][0][0]
    assert cmd[0:3] == ["diamond", "makedb", "--db"]

@pytest.mark.unit
def test_blastp(diamond, mock_successful_run, test_files):
    """Test blastp command."""
    # Test default format
    df = diamond.blastp(
        db=test_files['db'],
        query=test_files['fasta'],
        evalue=1e-10,
        threads=4
    )
    
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 1
    assert df.iloc[0]['qseqid'] == 'seq1'
    assert df.iloc[0]['pident'] == 98.5
    
    # Test custom format
    df = diamond.blastp(
        db=test_files['db'],
        query=test_files['fasta'],
        evalue=1e-10,
        threads=4,
        outfmt="6 qseqid sseqid pident"
    )
    
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 1
    assert list(df.columns) == ['qseqid', 'sseqid', 'pident']

@pytest.mark.unit
def test_blastx(diamond, mock_successful_run, test_files):
    """Test blastx command."""
    # Test default format
    df = diamond.blastx(
        db=test_files['db'],
        query=test_files['fasta'],
        evalue=1e-10,
        threads=4
    )
    
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 1
    assert df.iloc[0]['qseqid'] == 'seq1'
    assert df.iloc[0]['pident'] == 98.5
    
    # Test custom format
    df = diamond.blastx(
        db=test_files['db'],
        query=test_files['fasta'],
        evalue=1e-10,
        threads=4,
        outfmt="6 qseqid sseqid pident"
    )
    
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 1
    assert list(df.columns) == ['qseqid', 'sseqid', 'pident']

@pytest.mark.unit
def test_view(mocker, tmp_path):
    """Test view command."""
    mock_run = mocker.patch('subprocess.run')
    
    # Test tabular output with default format
    def mock_run_tabular(args, **kwargs):
        result = mocker.MagicMock()
        result.stdout = "seq1\tref1\t98.5\t100\t2\t0\t1\t100\t1\t100\t1e-50\t200.0\n"
        result.stderr = ""
        result.returncode = 0
        
        # Write to output file if specified
        for i, arg in enumerate(args):
            if arg == '--out' and i + 1 < len(args):
                with open(args[i + 1], 'w') as f:
                    f.write(result.stdout)
                break
        
        return result
    
    mock_run.side_effect = mock_run_tabular
    diamond = Diamond()
    
    result = diamond.view(
        daa="alignment.daa",
        outfmt=6
    )
    
    assert isinstance(result, pd.DataFrame)
    assert len(result) == 1
    assert result.iloc[0]['qseqid'] == 'seq1'
    assert result.iloc[0]['pident'] == 98.5
    
    # Test tabular output with custom format
    def mock_run_custom(args, **kwargs):
        result = mocker.MagicMock()
        result.stdout = "seq1\tref1\t98.5\n"
        result.stderr = ""
        result.returncode = 0
        
        # Write to output file if specified
        for i, arg in enumerate(args):
            if arg == '--out' and i + 1 < len(args):
                with open(args[i + 1], 'w') as f:
                    f.write(result.stdout)
                break
        
        return result
    
    mock_run.side_effect = mock_run_custom
    
    result = diamond.view(
        daa="alignment.daa",
        outfmt="6 qseqid sseqid pident"
    )
    
    assert isinstance(result, pd.DataFrame)
    assert len(result) == 1
    assert list(result.columns) == ['qseqid', 'sseqid', 'pident']
    
    # Test text output
    def mock_run_text(args, **kwargs):
        result = mocker.MagicMock()
        result.stdout = ">seq1\nACGT\n"
        result.stderr = ""
        result.returncode = 0
        
        # Write to output file if specified
        for i, arg in enumerate(args):
            if arg == '--out' and i + 1 < len(args):
                with open(args[i + 1], 'w') as f:
                    f.write(result.stdout)
                break
        
        return result
    
    mock_run.side_effect = mock_run_text
    
    # Text format (0)
    result = diamond.view(
        daa="alignment.daa",
        outfmt=0
    )
    
    assert isinstance(result, str)
    assert result == ">seq1\nACGT\n"

@pytest.mark.unit
def test_error_handling(diamond, mock_failed_run, test_files):
    """Test error handling."""
    with pytest.raises(DiamondError):
        diamond.makedb(
            db=test_files['db'],
            input_file=test_files['fasta']
        )

@pytest.mark.unit
def test_cluster(diamond, mock_successful_run, test_files):
    """Test cluster command."""
    df = diamond.cluster(
        db=test_files['db'],
        approx_id=90.0,
        threads=4
    )
    
    assert isinstance(df, pd.DataFrame)
    assert list(df.columns) == ['cluster_number', 'sequence_id', 'representative', 'cluster_size']

@pytest.mark.unit
def test_linclust(diamond, mock_successful_run, test_files):
    """Test linclust command."""
    df = diamond.linclust(
        db=test_files['db'],
        approx_id=90.0,
        threads=4
    )
    
    assert isinstance(df, pd.DataFrame)
    assert list(df.columns) == ['cluster_number', 'sequence_id', 'representative', 'cluster_size']

@pytest.mark.unit
def test_temp_file_cleanup(tmp_path):
    """Test temporary file cleanup."""
    diamond = Diamond()
    
    # Create a temporary file
    temp_file = diamond._create_temp_file()
    assert os.path.exists(temp_file)
    
    # Delete Diamond object
    del diamond
    
    # Check that temporary file was removed
    assert not os.path.exists(temp_file)

@pytest.mark.unit
def test_cluster_parser_mock_output():
    """Test ClusterParser.parse_clusters with a simple mock output."""
    from diamondonpy import ClusterParser
    
    # Create a simple mock output
    mock_output = "protein_1\tprotein_1\nprotein_2\tprotein_2\nprotein_3\tprotein_3\nprotein_4\tprotein_4\nprotein_4\tprotein_5"
    
    df = ClusterParser.parse_clusters(mock_output)
    
    # Verify that the result is a DataFrame with the expected columns
    expected_columns = {'cluster_number', 'sequence_id', 'representative', 'cluster_size'}
    assert isinstance(df, pd.DataFrame)
    assert set(df.columns) == expected_columns
    
    # There should be 5 unique sequence entries
    expected_ids = {"protein_1", "protein_2", "protein_3", "protein_4", "protein_5"}
    assert set(df['sequence_id']) == expected_ids
    
    # Group by cluster_number and check cluster sizes: one cluster should have 2 members and the rest 1
    cluster_sizes = df.groupby('cluster_number')['cluster_size'].first().tolist()
    assert sorted(cluster_sizes) == [1, 1, 1, 2]
    
    # For the cluster with size 2, verify that all representatives are the same
    for _, group in df.groupby('cluster_number'):
        if group['cluster_size'].iloc[0] == 2:
            rep = group['representative'].iloc[0]
            assert all(group['representative'] == rep)

@pytest.mark.unit
def test_validate_matrix_and_gaps(diamond):
    """Test matrix and gap penalty validation."""
    # Test with valid matrix and defaults
    options = {"matrix": "BLOSUM62"}
    result = diamond._validate_matrix_and_gaps(options)
    assert result["gapopen"] == 11
    assert result["gapextend"] == 1
    
    # Test with valid matrix and specified gap penalties
    options = {"matrix": "BLOSUM62", "gapopen": 9, "gapextend": 1}
    result = diamond._validate_matrix_and_gaps(options)
    assert result["gapopen"] == 9
    assert result["gapextend"] == 1
    
    # Test with invalid matrix
    with pytest.raises(ValueError, match="Invalid matrix"):
        diamond._validate_matrix_and_gaps({"matrix": "INVALID"})
    
    # Test with only one gap parameter
    with pytest.raises(ValueError, match="must also be specified"):
        diamond._validate_matrix_and_gaps({"matrix": "BLOSUM62", "gapopen": 9})
        
    # Test with invalid gap combination
    with pytest.raises(ValueError, match="Invalid gap penalties"):
        diamond._validate_matrix_and_gaps({"matrix": "BLOSUM62", "gapopen": 5, "gapextend": 5})

@pytest.mark.unit
def test_getseq(diamond, mock_successful_run, test_files):
    """Test getseq command."""
    # Test with output file
    output = diamond.getseq(
        db=test_files['db'],
        seq="1 2 3",
        out=test_files['out']
    )
    
    assert isinstance(output, str)
    assert ">seq1" in output
    assert "ACGT" in output
    
    # Test without output file (in-memory mode)
    diamond.use_memory = True
    output = diamond.getseq(
        db=test_files['db'],
        seq="1 2 3"
    )
    
    assert isinstance(output, str)
    assert ">seq1" in output
    assert "ACGT" in output
    
    # Reset memory mode
    diamond.use_memory = False

@pytest.mark.unit
def test_dbinfo(diamond, mock_successful_run, test_files):
    """Test dbinfo command."""
    output = diamond.dbinfo(db=test_files['db'])
    
    assert isinstance(output, str)
    assert ">seq1" in output
    assert "ACGT" in output

@pytest.mark.unit
def test_bidirectional_best_hit(diamond, mocker, test_files):
    """Test bidirectional best hit analysis."""
    # Mock blastp to return predictable DataFrames
    mock_blastp = mocker.patch.object(diamond, 'blastp')
    
    # Create mock results for both directions
    df1 = pd.DataFrame({
        'qseqid': ['seq1', 'seq1', 'seq2'],
        'sseqid': ['seq3', 'seq4', 'seq3'],
        'evalue': [0.001, 0.01, 0.0001]
    })
    
    df2 = pd.DataFrame({
        'qseqid': ['seq3', 'seq4', 'seq3'],
        'sseqid': ['seq1', 'seq1', 'seq2'],
        'evalue': [0.002, 0.02, 0.0002]
    })
    
    # Configure mock to return our test data
    mock_blastp.side_effect = [df1, df2]
    
    # Run the function
    result = diamond.bidirectional_best_hit(
        fasta1=test_files['db'],
        fasta2=test_files['db']
    )
    
    # Check the result is a DataFrame with expected structure
    assert isinstance(result, pd.DataFrame)
    assert 'qseqid' in result.columns
    assert 'sseqid' in result.columns
    
    # Check both directions of blastp were called
    assert mock_blastp.call_count == 2 

@pytest.mark.unit
def test_safe_command_string(diamond):
    """Test command string safety quoting."""
    cmd = ["diamond", "makedb", "--db", "test.dmnd", "--in", "file with spaces.fa"]
    safe_cmd = diamond._safe_command_string(cmd)
    
    # Ensure spaces in filenames are correctly quoted
    assert "file\\ with\\ spaces.fa" in safe_cmd or "'file with spaces.fa'" in safe_cmd 

@pytest.mark.unit
def test_prepare_fasta_input(diamond, test_files, tmp_path):
    """Test preparation of FASTA input."""
    # Test with existing file path
    input_path, use_stdin = diamond._prepare_fasta_input(test_files['fasta'])
    assert input_path == test_files['fasta']
    assert use_stdin is False
    
    # Test with string content in memory mode
    diamond.use_memory = True
    fasta_content = ">test\nACGT\n"
    input_path, use_stdin = diamond._prepare_fasta_input(fasta_content)
    assert input_path == fasta_content
    assert use_stdin is True
    
    # Test with string content in file mode
    diamond.use_memory = False
    input_path, use_stdin = diamond._prepare_fasta_input(fasta_content)
    assert os.path.exists(input_path)
    assert use_stdin is False
    assert input_path in diamond._temp_files
    
    # Reset memory mode
    diamond.use_memory = False 