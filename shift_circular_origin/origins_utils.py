"""
Utility functions for managing origin sequences in the origins.csv file.
"""

import csv
from pathlib import Path
from typing import List, Dict


def get_origins_file_path() -> Path:
    """
    Get the path to the origins.csv file.
    
    Returns:
        Path to origins.csv file
    """
    return Path(__file__).parent / "origins.csv"


def list_origins() -> List[Dict[str, str]]:
    """
    List all origin sequences from the origins.csv file.
    
    Returns:
        List of dictionaries with 'name' and 'sequence' keys
    """
    origins_file = get_origins_file_path()
    origins = []
    
    try:
        with open(origins_file, 'r', newline='', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            for row in reader:
                origins.append({
                    'name': row.get('name', '').strip(),
                    'sequence': row.get('sequence', '').strip()
                })
    except FileNotFoundError:
        print(f"Origins file not found: {origins_file}")
    except Exception as e:
        print(f"Error reading origins file: {e}")
    
    return origins


def add_origin(name: str, sequence: str) -> bool:
    """
    Add a new origin sequence to the origins.csv file.
    
    Args:
        name: Name of the origin
        sequence: DNA sequence
        
    Returns:
        True if successful, False otherwise
    """
    origins_file = get_origins_file_path()
    
    try:
        # Read existing origins
        existing_origins = list_origins()
        
        # Check if name already exists
        if any(origin['name'].lower() == name.lower() for origin in existing_origins):
            print(f"Origin '{name}' already exists. Use update_origin() to modify it.")
            return False
        
        # Add new origin
        existing_origins.append({
            'name': name.strip(),
            'sequence': sequence.strip().upper()
        })
        
        # Write back to file
        with open(origins_file, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=['name', 'sequence'])
            writer.writeheader()
            writer.writerows(existing_origins)
        
        print(f"Added origin '{name}' successfully.")
        return True
        
    except Exception as e:
        print(f"Error adding origin: {e}")
        return False


def remove_origin(name: str) -> bool:
    """
    Remove an origin sequence from the origins.csv file.
    
    Args:
        name: Name of the origin to remove
        
    Returns:
        True if successful, False otherwise
    """
    origins_file = get_origins_file_path()
    
    try:
        # Read existing origins
        existing_origins = list_origins()
        
        # Filter out the origin to remove
        filtered_origins = [
            origin for origin in existing_origins 
            if origin['name'].lower() != name.lower()
        ]
        
        if len(filtered_origins) == len(existing_origins):
            print(f"Origin '{name}' not found.")
            return False
        
        # Write back to file
        with open(origins_file, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=['name', 'sequence'])
            writer.writeheader()
            writer.writerows(filtered_origins)
        
        print(f"Removed origin '{name}' successfully.")
        return True
        
    except Exception as e:
        print(f"Error removing origin: {e}")
        return False


def main():
    """Main CLI entry point for origins management."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Manage origin sequences")
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # List command
    list_parser = subparsers.add_parser('list', help='List all origins')
    
    # Add command
    add_parser = subparsers.add_parser('add', help='Add new origin')
    add_parser.add_argument('name', help='Origin name')
    add_parser.add_argument('sequence', help='DNA sequence')
    
    # Remove command
    remove_parser = subparsers.add_parser('remove', help='Remove origin')
    remove_parser.add_argument('name', help='Origin name to remove')
    
    # Show file path command
    path_parser = subparsers.add_parser('path', help='Show origins.csv file path')
    
    args = parser.parse_args()
    
    if args.command == 'list':
        origins = list_origins()
        if origins:
            print(f"Found {len(origins)} origins:")
            for origin in origins:
                seq_preview = origin['sequence'][:50] + ('...' if len(origin['sequence']) > 50 else '')
                print(f"  {origin['name']}: {seq_preview}")
        else:
            print("No origins found.")
    
    elif args.command == 'add':
        add_origin(args.name, args.sequence)
    
    elif args.command == 'remove':
        remove_origin(args.name)
    
    elif args.command == 'path':
        print(f"Origins file: {get_origins_file_path()}")
    
    else:
        parser.print_help()


if __name__ == "__main__":
    main()