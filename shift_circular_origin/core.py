"""
Core functionality for rotating circular DNA sequences to start at specified origins.

This module contains the main SequenceRotator class and related functions for processing
DNA sequences to ensure they start at known origin sequences.
"""

import os
import re
from pathlib import Path
from typing import List, Tuple, Optional, Union
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class SequenceRotator:
    """
    A class for rotating circular DNA sequences to start at specified origin sequences.
    
    This class handles the rotation of circular DNA sequences (like plasmids) so that they
    begin at known origin sequences. It properly handles sequences where the origin spans
    the circular boundary (end-to-start wrap-around).
    """
    
    # Default candidate origin sequences (same as R script)
    DEFAULT_CANDIDATES = [
        "ATGAACAAGAGCGCCGCCGCTGGCCTGCTGGGCTATGCCCGCGTCAGCACCGACGACCAGGACTTGACCAACCAACGGGCCGAACTGCACGCGGCCGGCTGCACCAAGCTGTTTTCCGAGAAGATCACCGGCACCAGGCGCGACCGCCCG",  # pVS1 Long
        "CGTGCGGCTGCATGAAATCCTGGCCGGTTTGTCTGATGCCAAGCTGGCGGCCTGGCCGGCCAGCTTGGCCGCTGAAGAAACCGAGCGCCGCCGTCTAAAAAGGTGATGTGTATTTGAGTAAAACAGCTTGCGTCATGCGGTCGCTGCGTA",  # pVS1 short
        "TGGCGCCGGCCAGCGAGACGAGCAAGATTGGCCGCCGCCCGAAACGATCCGACAGCGCGCCCAGCACAGGTGCGCAGGCAAATTGCACCAACGCATACAGCGCCAGCAGAATGCCATAGTGGGCGGTGACGTCGTTCGAGTGAACCAGAT",  # pRNAiGG
        "CGCTGGCTGCTGAACCCCCAGCCGGAACTGACCCCACAAGGCCCTAGCGTTTGCAATGCACCAGGTCATCATTGACCCAGGCGTGTTCCACCAGGCCGCTGCCTCGCAACTCTTCGCAGGCTTCGCCGACCTGCTCGCGCCACTTCTTCA",  # RK2
        "GACGAAAGGGCCTCGTGATACGCCTATTTTTATAGGTTAATGTCATGATAATAATGGTTTCTTAGACGTCAGGTGGCACTTTTCGGGGAAATGTGCGCGGAACCCCTATTTGTTTATTTTTCTAAATACATTCAAATATGTATCCGCTCA",  # pL0V
        "AGGCACGAACCCAGTGGACATAAGCCTGTTCGGTTCGTAAGCTGTAATGCAAGTAGCGTATGCGCTCACGCAACTGGTCCAGAACCTTGACCGAACGCAGCGGTGGTAACGGCGCAGTGGCGGTTTTCATGGCTTGTTATGACTGTTTTT",  # pL0V uni
        "GGGTCCCCAATAATTACGATTTAAATTTGTGTCTCAAAATCTCTGATGTTACATTGCACAAGATAAAAATATATCATCATGAACAATAAAACTGTCTGCTTACATAAACAGTAATACAAGGGGTGTTATGAGCCATATTCAGCGTGAAAC",  # pSEVA23
        "GGGTCCCCAATAATTACGATTTAAATTTGACATAAGCCTGTTCGGTTCGTAAACTGTAATGCAAGTAGCGTATGCGCTCACGCAACTGGTCCAGAACCTTGACCGAACGCAGCGGTGGTAACGGCGCAGTGGCGGTTTTCATGGCTTGTT",  # pSEVA69
        "GCTTCACGCTGCCGCAAGCACTCAGGGCGCAAGGGCTGCTAAAGGAAGCGGAACACGTAGAAAGCCAGTCCGCAGAAACGGTGCTGACCCCGGATGAATGTCAGCTACTGGGCTATCTGGACAAGGGAAAACGCAAGCGCAAAGAGAAAG",  # pK18
        "ACGCGCCCTGTAGCGGCGCATTAAGCGCGGCGGGTGTGGTGGTTACGCGCAGCGTGACCGCTACACTTGCCAGCGCCCTAGCGCCCGCTCCTTTCGCTTTCTTCCCTTCCTTTCTCGCCACGTTCGCCGGCTTTCCCCGTCAAGCTCTAA",  # pET28
        "TCGCGTTATGCAGGCTTCCTCGCTCACTGACTCGCTGCGCTCGGTCGTTCGGCTGCGGCGAGCGGTATCAGCTCACTCAAAGGCGGTAATACGGTTATCCACAGAATCAGGGGATAACGCAGGAAAGAACATGTGAGCAAAAGGCCAGCA"   # pUAP1
    ]
    
    def __init__(self, candidate_sequences: Optional[List[str]] = None, max_length: int = 200000):
        """
        Initialize the SequenceRotator.
        
        Args:
            candidate_sequences: List of DNA sequences to search for as origins.
                                If None, uses default candidates.
            max_length: Maximum sequence length to process (default: 200,000 bp)
        """
        self.candidate_sequences = candidate_sequences or self.DEFAULT_CANDIDATES
        self.max_length = max_length
    
    def add_candidate_sequences(self, sequences: Union[str, List[str]]):
        """
        Add additional candidate sequences to search for.
        
        Args:
            sequences: Single sequence string or list of sequences to add.
                      Can also be a semicolon-separated string.
        """
        if isinstance(sequences, str):
            # Handle semicolon-separated sequences
            if ';' in sequences:
                new_sequences = [seq.strip().upper() for seq in sequences.split(';') if seq.strip()]
            else:
                new_sequences = [sequences.strip().upper()]
        else:
            new_sequences = [seq.strip().upper() for seq in sequences]
        
        # Add to existing candidates, avoiding duplicates
        for seq in new_sequences:
            if seq and seq not in self.candidate_sequences:
                self.candidate_sequences.append(seq)
    
    def _find_pattern_in_circular_sequence(self, sequence: str, pattern: str) -> Optional[int]:
        """
        Find a pattern in a circular sequence, handling boundary-spanning matches.
        
        This fixes the R script bug where matchPattern() treats sequences as linear,
        missing patterns that span the circular boundary (from end back to start).
        
        Args:
            sequence: DNA sequence string
            pattern: Pattern to search for
            
        Returns:
            Starting position of pattern (0-based) or None if not found
        """
        # Create doubled sequence to handle circular boundary cases
        # Example: "AAATTTCCCGGG" becomes "AAATTTCCCGGGAAATTTCCCGGG"
        # This allows "GGGAAA" to be found as a contiguous substring
        doubled_sequence = sequence + sequence
        
        # Search in the region that could contain the pattern
        # We only need to search up to len(sequence) + len(pattern) - 1
        # to catch all possible boundary-spanning matches
        search_region = doubled_sequence[:len(sequence) + len(pattern) - 1]
        
        match = re.search(pattern, search_region)
        if match:
            return match.start()
        
        return None
    
    def _rotate_sequence(self, sequence: str, start_pos: int) -> str:
        """
        Rotate a circular sequence to start at the specified position.
        
        Args:
            sequence: Original DNA sequence
            start_pos: Position to start the rotated sequence (0-based)
            
        Returns:
            Rotated sequence string
        """
        if start_pos == 0:
            return sequence
        
        return sequence[start_pos:] + sequence[:start_pos]
    
    def _reverse_complement(self, sequence: str) -> str:
        """
        Generate the reverse complement of a DNA sequence.
        
        Args:
            sequence: DNA sequence string
            
        Returns:
            Reverse complement sequence string
        """
        complement_map = str.maketrans('ATCG', 'TAGC')
        return sequence.translate(complement_map)[::-1]
    
    def process_sequence(self, seq_record: SeqRecord) -> SeqRecord:
        """
        Process a single sequence record to rotate it to start at a candidate origin.
        
        This method implements the core logic:
        1. Check if sequence already starts with a candidate
        2. Search forward strand for candidates
        3. If not found, search reverse complement
        4. Rotate sequence to start at found candidate
        
        Args:
            seq_record: BioPython SeqRecord object
            
        Returns:
            Processed SeqRecord with potentially rotated sequence
        """
        sequence = str(seq_record.seq).upper()
        
        # Skip sequences that are too long
        if len(sequence) > self.max_length:
            return seq_record
        
        # Check if sequence already starts with a candidate
        for candidate in self.candidate_sequences:
            if sequence.startswith(candidate):
                return seq_record
        
        # Search forward strand
        for candidate in self.candidate_sequences:
            start_pos = self._find_pattern_in_circular_sequence(sequence, candidate)
            if start_pos is not None:
                if start_pos == 0:
                    return seq_record
                else:
                    rotated_seq = self._rotate_sequence(sequence, start_pos)
                    new_record = SeqRecord(
                        Seq(rotated_seq),
                        id=seq_record.id,
                        description=seq_record.description
                    )
                    return new_record
        
        # Search reverse complement
        rev_comp_sequence = self._reverse_complement(sequence)
        for candidate in self.candidate_sequences:
            start_pos = self._find_pattern_in_circular_sequence(rev_comp_sequence, candidate)
            if start_pos is not None:
                if start_pos == 0:
                    new_record = SeqRecord(
                        Seq(rev_comp_sequence),
                        id=seq_record.id,
                        description=seq_record.description
                    )
                    return new_record
                else:
                    rotated_seq = self._rotate_sequence(rev_comp_sequence, start_pos)
                    new_record = SeqRecord(
                        Seq(rotated_seq),
                        id=seq_record.id,
                        description=seq_record.description
                    )
                    return new_record
        
        # No candidate found, return unchanged
        return seq_record
    
    def process_fasta_file(self, input_file: Union[str, Path]) -> List[SeqRecord]:
        """
        Process all sequences in a FASTA file.
        
        Args:
            input_file: Path to input FASTA file
            
        Returns:
            List of processed SeqRecord objects
        """
        processed_records = []
        
        try:
            for record in SeqIO.parse(input_file, "fasta"):
                processed_record = self.process_sequence(record)
                processed_records.append(processed_record)
        except Exception as e:
            print(f"Error processing file {input_file}: {e}")
            return []
        
        return processed_records
    
    def process_folder(self, folder_path: Union[str, Path], output_file: Optional[Union[str, Path]] = None) -> bool:
        """
        Process all FASTA files in a folder and save results.
        
        Args:
            folder_path: Path to folder containing FASTA files
            output_file: Path for output file. If None, saves as 'origin_specified_sequences.fasta' in input folder
            
        Returns:
            True if successful, False otherwise
        """
        folder_path = Path(folder_path)
        if not folder_path.exists():
            print(f"Error: Folder {folder_path} does not exist")
            return False
        
        # Find FASTA files
        fasta_patterns = ["*.fa", "*.fasta", "*.fas"]
        fasta_files = []
        for pattern in fasta_patterns:
            fasta_files.extend(folder_path.glob(pattern))
        
        if not fasta_files:
            print(f"No FASTA files found in {folder_path}")
            return False
        
        print(f"Found {len(fasta_files)} FASTA files")
        
        # Process all files
        all_processed_records = []
        for fasta_file in fasta_files:
            print(f"Processing: {fasta_file.name}")
            processed_records = self.process_fasta_file(fasta_file)
            all_processed_records.extend(processed_records)
        
        if not all_processed_records:
            print("No sequences were processed")
            return False
        
        # Set output file path
        if output_file is None:
            output_file = folder_path / "origin_specified_sequences.fasta"
        else:
            output_file = Path(output_file)
        
        # Write output
        try:
            SeqIO.write(all_processed_records, output_file, "fasta")
            print(f"Output written to: {output_file}")
            print(f"Processed {len(all_processed_records)} sequences")
            return True
        except Exception as e:
            print(f"Error writing output file: {e}")
            return False


def rotate_sequences_in_folder(folder_path: Union[str, Path], 
                              additional_candidates: Optional[Union[str, List[str]]] = None,
                              output_file: Optional[Union[str, Path]] = None,
                              max_length: int = 200000) -> bool:
    """
    Convenience function to process all FASTA files in a folder.
    
    Args:
        folder_path: Path to folder containing FASTA files
        additional_candidates: Additional origin sequences to search for (semicolon-separated string or list)
        output_file: Path for output file (optional)
        max_length: Maximum sequence length to process (default: 200,000 bp)
        
    Returns:
        True if successful, False otherwise
    
    Example:
        # Basic usage
        rotate_sequences_in_folder("/path/to/fasta/files")
        
        # With additional origin sequences
        rotate_sequences_in_folder("/path/to/fasta/files", 
                                  additional_candidates="ATGCCCGGG;CCCGGGATG")
    """
    rotator = SequenceRotator(max_length=max_length)
    
    if additional_candidates:
        rotator.add_candidate_sequences(additional_candidates)
    
    return rotator.process_folder(folder_path, output_file)