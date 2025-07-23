# Shift Circular Origin

A Python package for rotating circular DNA sequences to start at specified origin sequences. This tool is particularly useful for standardizing plasmids and other circular genetic elements by ensuring they begin at known origin of replication sequences.

## Features

- **Circular Boundary Handling**: Properly detects origin sequences that span the circular boundary (wrap from end to start)
- **Dual Strand Search**: Searches both forward and reverse complement strands
- **Multiple Origin Support**: Includes built-in database of common plasmid origins and supports custom sequences
- **Batch Processing**: Process entire folders of FASTA files at once
- **Length Filtering**: Configurable maximum sequence length (default: 200kb)
- **Flexible Input**: Supports semicolon-separated custom origin sequences

## Installation

### From PyPI (when published)
```bash
pip install shift-circular-origin
```

### From Source
```bash
git clone https://github.com/yourusername/shift-circular-origin.git
cd shift-circular-origin
pip install -e .
```

## Quick Start

### Command Line Usage

Process all FASTA files in a folder:
```bash
shift-circular-origin /path/to/fasta/files
```

Add custom origin sequences:
```bash
shift-circular-origin /path/to/fasta/files --additional-origins "ATGCCCGGG;CCCGGGATG"
```

Specify output file:
```bash
shift-circular-origin /path/to/fasta/files --output results.fasta
```

### Python API Usage

```python
from shift_circular_origin import SequenceRotator, rotate_sequences_in_folder

# Simple batch processing
rotate_sequences_in_folder("/path/to/fasta/files")

# With additional origin sequences
rotate_sequences_in_folder(
    "/path/to/fasta/files", 
    additional_candidates="ATGCCCGGG;CCCGGGATG"
)

# Advanced usage with custom rotator
rotator = SequenceRotator()
rotator.add_candidate_sequences(["ATGCCCGGG", "CCCGGGATG"])
rotator.process_folder("/path/to/fasta/files", "output.fasta")
```

## How It Works

1. **Input Processing**: Reads all FASTA files from the specified folder
2. **Length Filtering**: Excludes sequences longer than the specified maximum (default: 200kb)
3. **Origin Detection**: Searches for candidate origin sequences in both forward and reverse complement orientations
4. **Circular Boundary Handling**: Uses a doubled-sequence approach to detect origins that span the end-to-start boundary
5. **Sequence Rotation**: Rotates matching sequences so they start at the detected origin
6. **Output Generation**: Writes all processed sequences to a single FASTA file

## Built-in Origin Sequences

The package includes origins from common plasmid backbones:

- **pVS1**: Long and short variants
- **pRNAiGG**: Gateway-compatible vectors
- **RK2**: Broad-host-range plasmids
- **pL0V**: SEVA collection vectors
- **pSEVA23/69**: Standard biological parts
- **pK18**: Cloning vectors
- **pET28**: Expression vectors
- **pUAP1**: Plant transformation vectors

## Technical Details

### Circular Boundary Problem

The original R implementation had a critical bug where origin sequences spanning the circular boundary (from end back to start) were not detected. This Python implementation fixes this by:

1. Creating a doubled sequence (original + original)
2. Searching for patterns in the first copy plus pattern length
3. Properly handling matches that extend beyond the original sequence length

### Algorithm Steps

1. Check if sequence already starts with any candidate origin
2. Search forward strand using circular-aware pattern matching
3. If found, rotate sequence to start at that position
4. If not found, search reverse complement strand
5. If found in reverse complement, use the reverse complement and rotate
6. If no origin found, leave sequence unchanged

## API Reference

### SequenceRotator Class

Main class for handling sequence rotation operations.

**Constructor:**
```python
SequenceRotator(candidate_sequences=None, max_length=200000)
```

**Key Methods:**
- `add_candidate_sequences(sequences)`: Add custom origin sequences
- `process_sequence(seq_record)`: Process a single BioPython SeqRecord
- `process_fasta_file(input_file)`: Process all sequences in a FASTA file
- `process_folder(folder_path, output_file)`: Batch process folder of FASTA files

### Convenience Function

```python
rotate_sequences_in_folder(folder_path, additional_candidates=None, output_file=None, max_length=200000)
```

## Requirements

- Python >= 3.7
- BioPython >= 1.70

## Project Structure

```
shift-circular-origin/
├── shift_circular_origin/
│   ├── __init__.py          # Package initialization
│   ├── core.py              # Main functionality
│   └── cli.py               # Command-line interface
├── setup.py                 # Package configuration
├── requirements.txt         # Dependencies
├── README.md               # This file
├── LICENSE                 # MIT License
└── .gitignore              # Git ignore rules
```

## Contributing

1. Fork the repository
2. Create a feature branch: `git checkout -b feature-name`
3. Make your changes and add tests
4. Run tests: `python -m pytest`
5. Commit changes: `git commit -am 'Add feature'`
6. Push to branch: `git push origin feature-name`
7. Submit a pull request

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use this software in your research, please cite:

```
Shift Circular Origin: A Python package for rotating circular DNA sequences
URL: https://github.com/JKourelis/shift-circular-origin
```

## Changelog

### Version 1.0.0
- Initial release
- Fixed circular boundary detection bug from original R implementation
- Added command-line interface
- Added support for custom origin sequences
- Comprehensive error handling and logging