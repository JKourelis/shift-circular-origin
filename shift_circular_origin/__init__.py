"""
Shift Circular Origin - DNA Sequence Origin Rotation Package

This package provides functionality to rotate circular DNA sequences so that they start
at specified origin sequences. It handles plasmids and other circular genetic elements
by searching for known origin sequences and rotating the sequence to begin at that point.

Key Features:
- Handles circular sequence boundaries (sequences spanning start/end)
- Searches both forward and reverse complement strands
- Supports multiple candidate origin sequences
- Filters sequences by length (default: <= 200kb)
- Batch processing of FASTA files

Classes:
    SequenceRotator: Main class for handling sequence rotation operations

Functions:
    rotate_sequences_in_folder: Convenience function for batch processing
"""

from .core import SequenceRotator, rotate_sequences_in_folder

__version__ = "1.0.0"
__author__ = "Your Name"
__email__ = "your.email@example.com"

__all__ = ["SequenceRotator", "rotate_sequences_in_folder"]