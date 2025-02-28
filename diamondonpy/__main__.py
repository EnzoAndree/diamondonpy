import sys
import subprocess
from .diamond import Diamond
from .errors import DiamondError

def main():
    """
    Command line interface for diamondonpy.
    
    Usage: diamondonpy <subcommand> [options...]
    """
    if len(sys.argv) < 2:
        print("Usage: diamondonpy <subcommand> [options...]")
        sys.exit(1)
        
    subcommand = sys.argv[1]
    diamond = Diamond()
    
    try:
        # Pass through all arguments to the diamond executable
        command = [diamond.executable, subcommand] + sys.argv[2:]
        result = subprocess.run(command, capture_output=True, text=True, check=True)
        print(result.stdout)
    except subprocess.CalledProcessError as e:
        print(f"Error: {e.stderr}", file=sys.stderr)
        sys.exit(1)
    except DiamondError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main() 