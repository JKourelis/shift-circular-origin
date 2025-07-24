"""
Example usage of the shift-circular-origin package.

This script demonstrates how to use the package programmatically.
"""

from shift_circular_origin import SequenceRotator, rotate_sequences_in_folder
from shift_circular_origin.origins_utils import list_origins, add_origin, get_origins_file_path
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def example_basic_usage():
    """
    Example of basic usage with the convenience function.
    """
    print("=== Basic Usage Example ===")
    
    # Simple usage - process all FASTA files in a folder
    folder_path = "/path/to/your/fasta/files"  # Change this path
    success = rotate_sequences_in_folder(folder_path)
    
    if success:
        print("Successfully processed sequences!")
    else:
        print("Processing failed - check folder path and file contents")

def example_origins_management():
    """
    Example of managing origin sequences via CSV.
    """
    print("=== Origins Management Example ===")
    
    # Show current origins
    print("Current origins:")
    origins = list_origins()
    for origin in origins[:3]:  # Show first 3
        seq_preview = origin['sequence'][:30] + "..."
        print(f"  {origin['name']}: {seq_preview}")
    
    print(f"\nTotal origins loaded: {len(origins)}")
    
    # Show where the CSV file is located
    csv_path = get_origins_file_path()
    print(f"Origins CSV file: {csv_path}")
    
    # Example of adding a custom origin (commented out to avoid modifying the file)
    # add_origin("my_custom_origin", "ATGCCCGGGAAATTTGCGATCGC")
    # print("Added custom origin!")

def example_with_custom_origins():
    """
    Example using custom origin sequences via command line.
    """
    print("=== Custom Origins Example ===")
    
    # Add custom origin sequences (semicolon-separated)
    custom_origins = "ATGCCCGGGAAATTT;CCCGGGATGAAATTT;GGGAAATTTCCCGGG"
    folder_path = "/path/to/your/fasta/files"  # Change this path
    output_file = "custom_origin_results.fasta"
    
    success = rotate_sequences_in_folder(
        folder_path,
        additional_candidates=custom_origins,
        output_file=output_file
    )
    
    if success:
        print(f"Results saved to {output_file}")

def example_circular_boundary_fix():
    """
    Example demonstrating the circular boundary fix.
    
    This shows how the package handles origin sequences that span
    the circular boundary (end-to-start wrap).
    """
    print("=== Circular Boundary Fix Example ===")
    
    # Create a sequence where the origin spans the boundary
    # Let's say our origin is "ATGCCC" but the sequence is stored as:
    # "CCCGGGAAATTTGCGATCGCTAGCTATG"
    # The origin "ATGCCC" spans from position 24 ("ATG") to position 0-2 ("CCC")
    
    boundary_spanning_seq = "CCCGGGAAATTTGCGATCGCTAGCTATG"
    origin_to_find = "ATGCCC"
    
    rotator = SequenceRotator()
    rotator.add_candidate_sequences([origin_to_find])
    
    seq_record = SeqRecord(
        Seq(boundary_spanning_seq), 
        id="boundary_test", 
        description="Sequence with boundary-spanning origin"
    )
    
    processed_record = rotator.process_sequence(seq_record)
    
    print(f"Original sequence: {boundary_spanning_seq}")
    print(f"Looking for origin: {origin_to_find}")
    print(f"Processed sequence: {str(processed_record.seq)}")
    print(f"Starts with origin: {str(processed_record.seq).startswith(origin_to_find)}")

if __name__ == "__main__":
    print("Shift Circular Origin - Usage Examples")
    print("=====================================")
    
    # Run examples (comment out the ones that require actual files)
    # example_basic_usage()
    # example_with_custom_origins()
    
    # These work without external files
    example_origins_management()
    example_circular_boundary_fix()
    
    print("\nFor the file-based examples, update the folder_path variables with actual paths.")
    print("Install the package first: pip install -e .")
    print("\nTo manage origins:")
    print("  origins-manager list")
    print("  origins-manager add my_origin ATGCCCGGG")