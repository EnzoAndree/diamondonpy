"""Pytest configuration and fixtures."""
import pytest
from pathlib import Path
import subprocess
from diamondonpy import Diamond

@pytest.fixture
def mock_successful_run(mocker):
    """Mock a successful subprocess run."""
    mock = mocker.patch('subprocess.run')
    
    def mock_run(args, **kwargs):
        result = mocker.MagicMock()
        result.stderr = ""
        result.returncode = 0
        
        # Handle version check
        if args[1] == "--version":
            result.stdout = "diamond version 2.1.0"
            return result
            
        # Handle BLAST-like outputs with output file
        if any(cmd in args for cmd in ['blastp', 'blastx', 'view']):
            # Get output format if specified
            outfmt = None
            for i, arg in enumerate(args):
                if arg == '--outfmt' and i + 1 < len(args):
                    outfmt = args[i + 1]
                    break
            
            # Get output file if specified
            out_file = None
            for i, arg in enumerate(args):
                if arg == '--out' and i + 1 < len(args):
                    out_file = args[i + 1]
                    break
            
            # Generate mock data based on full column set
            mock_data = {
                'qseqid': 'seq1',
                'sseqid': 'ref1',
                'pident': '98.5',
                'length': '100',
                'mismatch': '2',
                'gapopen': '0',
                'qstart': '1',
                'qend': '100',
                'sstart': '1',
                'send': '100',
                'evalue': '1e-50',
                'bitscore': '200.0'
            }
            
            # Format output according to outfmt
            if isinstance(outfmt, str) and outfmt.startswith('6 '):
                columns = outfmt.split(' ')[1:]
                if columns:  # Only if we have columns specified
                    output_parts = []
                    for col in columns:
                        output_parts.append(mock_data.get(col, ''))
                    output = '\t'.join(output_parts) + '\n'
                else:
                    # In case of empty column specification
                    output = "seq1\tref1\t98.5\t100\t2\t0\t1\t100\t1\t100\t1e-50\t200.0\n"
            else:
                # Default format (6 or other)
                output = "seq1\tref1\t98.5\t100\t2\t0\t1\t100\t1\t100\t1e-50\t200.0\n"
            
            if out_file:
                with open(out_file, 'w') as f:
                    f.write(output)
            result.stdout = output
            return result
            
        # Handle clustering outputs
        if any(cmd in args for cmd in ['cluster', 'linclust']):
            result.stdout = "seq1\tseq2\t98.5\n"
            return result

        # Handle getseq and dbinfo output
        if any(cmd in args for cmd in ['getseq', 'dbinfo']):
            # Mock output for getseq
            output = ">seq1\nACGT\n>seq2\nTGCA\n"
            result.stdout = output
            
            # Check if output file is specified
            out_file = None
            for i, arg in enumerate(args):
                if arg == '--out' and i + 1 < len(args):
                    out_file = args[i + 1]
                    break
                    
            if out_file:
                with open(out_file, 'w') as f:
                    f.write(output)
            
            return result
            
        # Default output
        result.stdout = "mock output"
        return result
        
    mock.side_effect = mock_run
    return mock

@pytest.fixture
def mock_failed_run(mocker):
    """Mock a failed subprocess run."""
    mock = mocker.patch('subprocess.run')
    
    def mock_run(args, **kwargs):
        if args[1] == "--version":
            raise subprocess.CalledProcessError(1, args, output="", stderr="mock error")
        if kwargs.get('check', False):
            raise subprocess.CalledProcessError(1, args, output="", stderr="mock error")
        result = mocker.MagicMock()
        result.stdout = ""
        result.stderr = "mock error"
        result.returncode = 1
        return result
        
    mock.side_effect = mock_run
    return mock

@pytest.fixture
def diamond(mock_successful_run):
    """Return a Diamond instance for testing."""
    return Diamond()

@pytest.fixture
def test_files(tmp_path):
    """Create temporary test files."""
    # Create a mock FASTA file
    fasta_file = tmp_path / "test.fasta"
    fasta_file.write_text(">seq1\nMKLVVVPET\n>seq2\nMIAVAVALA\n")
    
    # Create output paths
    db_file = tmp_path / "test.dmnd"
    out_file = tmp_path / "results.txt"
    
    return {
        'fasta': str(fasta_file),
        'db': str(db_file),
        'out': str(out_file)
    }