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

```bash
git clone https://github.com/yourusername/shift-circular-origin.git
cd shift-circular-origin
pip install -e .
```

## Quick Start

### Command Line Usage

Process all FASTA files in a folder:
```bash
shift-circular-origin /path/to/fasta/folder
```

Add custom origin sequences:
```bash
shift-circular-origin /path/to/fasta/folder --additional-origins "ATGCCCGGG;CCCGGGATG"
```

Specify output file:
```bash
shift-circular-origin /path/to/fasta/folder --output results.fasta
```

### Python Usage

```python
from shift_circular_origin import rotate_sequences_in_folder

# Basic usage
rotate_sequences_in_folder("fasta_folder")

# With custom origins
rotate_sequences_in_folder("fasta_folder", additional_candidates="ATGCCCGGG;CCCGGGATG")
```

## How It Works

1. **Input Processing**: Reads all FASTA files from the specified folder
2. **Length Filtering**: Excludes sequences longer than the specified maximum (default: 200kb)
3. **Origin Detection**: Searches for candidate origin sequences in both forward and reverse complement orientations
4. **Circular Sequence Handling**: Uses a doubled-sequence approach to detect origins that span the end-to-start boundary
5. **Sequence Rotation**: Rotates matching sequences so they start at the detected origin
6. **Output Generation**: Writes all processed sequences to a single FASTA file

## Built-in Origin Sequences

The package includes origins from common plasmid backbones, defined in `origins.csv`:

- **pVS1**: Long and short variants
- **pRNAiGG**: Gateway-compatible vectors
- **RK2**: Broad-host-range plasmids
- **pL0V**: SEVA collection vectors
- **pSEVA23/69**: Standard biological parts
- **pK18**: Cloning vectors
- **pET28**: Expression vectors
- **pUAP1**: Plant transformation vectors

## Customizing Origin Sequences

You can customize the origin sequences by editing the `origins.csv` file in the installed package.

### Using the Origins Manager CLI

```bash
# List all current origins
origins-manager list

# Find the origins.csv file location
origins-manager path

# Add a new origin
origins-manager add my_origin ATGCCCGGGAAATTT

# Remove an origin
origins-manager remove my_origin
```

### Manual Editing

1. **Find the file location:**
```python
import shift_circular_origin
from pathlib import Path
csv_path = Path(shift_circular_origin.__file__).parent / "origins.csv"
print(csv_path)
```

2. **Edit the CSV file:**
```csv
name,sequence
my_custom_origin,ATGCCCGGGAAATTT
pVS1_long,ATGAACAAGAGCGCC...
```

3. **Restart Python** to load the new origins

The CSV format is simple: `name,sequence` with one origin per line.

## Requirements

- Python >= 3.7
- BioPython >= 1.70

## Project Structure

```
shift-circular-origin/
├── shift_circular_origin/
│   ├── __init__.py          # Package initialization
│   ├── core.py              # Main functionality
│   ├── cli.py               # Command-line interface
│   ├── origins_utils.py     # Origins management utilities
│   └── origins.csv          # Origin sequences database (editable)
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
URL: https://github.com/yourusername/shift-circular-origin
```

## Changelog

### Version 1.0.0
- Initial release
- Fixed circular boundary detection bug from original R implementation
- Added command-line interface
- Added support for custom origin sequences
- Comprehensive error handling and logging