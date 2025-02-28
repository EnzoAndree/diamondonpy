"""Unit tests for the command line interface."""
import sys
import pytest
import subprocess
from diamondonpy.__main__ import main
from diamondonpy.errors import DiamondError

def test_main_no_args(capsys):
    """Test main function with no arguments."""
    with pytest.raises(SystemExit) as exc_info:
        sys.argv = ["diamondonpy"]
        main()
    
    assert exc_info.value.code == 1
    captured = capsys.readouterr()
    assert "Usage: diamondonpy <subcommand> [options...]" in captured.out

def test_main_with_args(mocker, capsys):
    """Test main function with arguments."""
    # Mock subprocess.run for both version check and command
    mock_run = mocker.patch('subprocess.run')
    
    def mock_run_func(args, **kwargs):
        result = mocker.MagicMock()
        result.returncode = 0
        if args[1] == "--version":
            result.stdout = "diamond version 2.1.0"
        else:
            result.stdout = "mock output"
        result.stderr = ""
        return result
    
    mock_run.side_effect = mock_run_func
    
    sys.argv = ["diamondonpy", "makedb", "--db", "test.dmnd", "--in", "test.fasta"]
    main()
    
    # Verify both calls were made
    assert mock_run.call_count == 2
    assert mock_run.call_args_list[0][0][0] == ["diamond", "--version"]
    assert mock_run.call_args_list[1][0][0] == ["diamond", "makedb", "--db", "test.dmnd", "--in", "test.fasta"]
    
    captured = capsys.readouterr()
    assert "mock output" in captured.out

def test_main_error(mocker, capsys):
    """Test main function with error."""
    # Mock subprocess.run to fail on version check
    mock_run = mocker.patch('subprocess.run')
    
    def mock_run_func(args, **kwargs):
        if args[1] == "--version":
            result = mocker.MagicMock()
            result.stdout = "diamond version 2.1.0"
            result.stderr = ""
            result.returncode = 0
            return result
        raise subprocess.CalledProcessError(1, args, output="", stderr="mock error")
    
    mock_run.side_effect = mock_run_func
    
    with pytest.raises(SystemExit) as exc_info:
        sys.argv = ["diamondonpy", "makedb", "--db", "test.dmnd"]
        main()
    
    assert exc_info.value.code == 1
    captured = capsys.readouterr()
    assert "Error: mock error" in captured.err 

@pytest.mark.unit
def test_main_unexpected_error(mocker, capsys):
    """Test main function when an unexpected exception occurs."""
    # Import main function
    from diamondonpy.__main__ import main
    
    # Mock Diamond to return a mock object that doesn't fail
    mock_diamond = mocker.patch('diamondonpy.__main__.Diamond')
    diamond_instance = mocker.MagicMock()
    mock_diamond.return_value = diamond_instance
    
    # Mock sys.argv
    mocker.patch.object(sys, 'argv', ['diamondonpy', 'blastp'])
    
    # Mock subprocess.run to raise a generic exception
    mock_run = mocker.patch('subprocess.run')
    mock_run.side_effect = Exception("Unexpected test error")
    
    # Call main
    with pytest.raises(SystemExit) as e:
        main()
    
    # Check exit code
    assert e.value.code == 1
    
    # Check error message
    captured = capsys.readouterr()
    assert "Unexpected error: Unexpected test error" in captured.err

@pytest.mark.unit
def test_main_diamond_error(mocker, capsys):
    """Test main function when a DiamondError occurs."""
    # Import main function and error
    from diamondonpy.__main__ import main
    from diamondonpy.errors import DiamondError
    
    # Mock Diamond to return a mock object that doesn't fail
    mock_diamond = mocker.patch('diamondonpy.__main__.Diamond')
    diamond_instance = mocker.MagicMock()
    mock_diamond.return_value = diamond_instance
    
    # Mock sys.argv
    mocker.patch.object(sys, 'argv', ['diamondonpy', 'blastp'])
    
    # Mock subprocess.run to raise a DiamondError
    mock_run = mocker.patch('subprocess.run')
    mock_run.side_effect = DiamondError("Diamond test error")
    
    # Call main
    with pytest.raises(SystemExit) as e:
        main()
    
    # Check exit code
    assert e.value.code == 1
    
    # Check error message
    captured = capsys.readouterr()
    assert "Error: Diamond test error" in captured.err 