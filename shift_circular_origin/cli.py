"""
Command-line interface for shift-circular-origin package.
"""

import argparse
import sys
from pathlib import Path
from . import rotate_sequences_in_folder


def main():
    """
    Main CLI entry point for the shift-circular-origin package.
    """
    parser = argparse.ArgumentParser(
        description="Rotate circular DNA sequences to start at specified origin sequences",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process all FASTA files in a folder
  shift-circular-origin /path/to/fasta/files

  # Add custom origin sequences (semicolon-separated)
  shift-circular-origin /path/to/fasta/files --additional-origins "ATGCCCGGG;CCCGGGATG"

  # Specify output file
  shift-circular-origin /path/to/fasta/files --output results.fasta

  # Change maximum sequence length filter
  shift-circular-origin /path/to/fasta/files --max-length 500000
        """
    )
    
    parser.add_argument(
        "folder_path",
        help="Path to folder containing FASTA files to process"
    )
    
    parser.add_argument(
        "--additional-origins", "-a",
        help="Additional origin sequences to search for (semicolon-separated)",
        default=None
    )
    
    parser.add_argument(
        "--output", "-o",
        help="Output FASTA file path (default: origin_specified_sequences.fasta in input folder)",
        default=None
    )
    
    parser.add_argument(
        "--max-length", "-m",
        type=int,
        help="Maximum sequence length to process in base pairs (default: 200000)",
        default=200000
    )
    
    parser.add_argument(
        "--version",
        action="version",
        version="shift-circular-origin 1.0.0"
    )
    
    args = parser.parse_args()
    
    # Validate input folder
    folder_path = Path(args.folder_path)
    if not folder_path.exists():
        print(f"Error: Folder '{folder_path}' does not exist", file=sys.stderr)
        sys.exit(1)
    
    if not folder_path.is_dir():
        print(f"Error: '{folder_path}' is not a directory", file=sys.stderr)
        sys.exit(1)
    
    # Process sequences
    print(f"Processing FASTA files in: {folder_path}")
    if args.additional_origins:
        print(f"Additional origin sequences: {args.additional_origins}")
    if args.output:
        print(f"Output file: {args.output}")
    
    success = rotate_sequences_in_folder(
        folder_path=folder_path,
        additional_candidates=args.additional_origins,
        output_file=args.output,
        max_length=args.max_length
    )
    
    if success:
        print("Processing completed successfully!")
        sys.exit(0)
    else:
        print("Processing failed!", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()