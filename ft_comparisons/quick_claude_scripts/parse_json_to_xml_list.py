"""
Module for extracting vasprun.xml file paths from training/validation JSON specifications.

This module handles parsing of JSON files that specify DFT calculation data locations
and generates complete file paths for all relevant vasprun.xml files.
"""

import json
import re
from pathlib import Path
from typing import List, Dict, Any, Optional


def normalize_path(file_path: str, path_type: str = "auto") -> str:
    """
    Normalize file paths by extracting only the relevant portion after xml_collection or xml_set.
    
    Args:
        file_path: Original file path from JSON
        path_type: Type of path ("equil", "trainval", or "auto" for auto-detection)
    
    Returns:
        Normalized path string
    """
    if not file_path:
        return ""
    
    # Check for xml_collection
    if "xml_collection" in file_path:
        parts = file_path.split("xml_collection/")
        return parts[1]
    
    # Check for xml_set
    elif "xml_set" in file_path:
        parts = file_path.split("xml_set/")
        return parts[1]
    
    # If neither marker is found, return original path
    return file_path


def parse_file_nums(file_nums_str: str) -> List[int]:
    """
    Parse file number specifications into a list of integers.
    
    Args:
        file_nums_str: String like "1..51" or "*"
    
    Returns:
        List of file numbers (empty list for "*" which means all files)
    """
    if file_nums_str == "*":
        return []  # Empty list indicates all files
    
    # Handle range format like "1..51"
    match = re.match(r"(\d+)\.\.(\d+)", file_nums_str)
    if match:
        start = int(match.group(1))
        end = int(match.group(2))
        return list(range(start, end + 1))
    
    # Handle single number
    if file_nums_str.isdigit():
        return [int(file_nums_str)]
    
    raise ValueError(f"Invalid file_nums format: {file_nums_str}")


def extract_equil_paths(system: Dict[str, Any]) -> List[str]:
    """
    Extract equilibrium file paths from a system specification.
    
    Args:
        system: System dictionary containing equil specifications
    
    Returns:
        List of normalized file paths
    """
    paths = []
    
    for equil_entry in system.get("equil", []):
        xml_file = equil_entry.get("xml_file", "")
        if xml_file:  # Skip empty strings
            normalized = normalize_path(xml_file, "equil")
            paths.append(normalized)
    
    return paths


def extract_trainval_paths(trainval_entry: Dict[str, Any], base_dir: Optional[str] = None) -> List[str]:
    """
    Extract training/validation file paths from a single trainval entry.
    
    Args:
        trainval_entry: Dictionary containing xml_dir and file_nums
        base_dir: Base directory to prepend to normalized paths for globbing
    
    Returns:
        List of normalized file paths
    """
    paths = []
    xml_dir = trainval_entry.get("xml_dir", "")
    file_nums_list = trainval_entry.get("file_nums", [])
    
    if not xml_dir:
        return paths
    
    # Normalize the directory path
    normalized_dir = normalize_path(xml_dir, "trainval")
    
    # Handle each file_nums specification
    for file_nums_spec in file_nums_list:
        if file_nums_spec == "*":
            # For wildcard, glob the directory to find all vasprun*.xml files
            if base_dir:
                search_path = Path(base_dir) / normalized_dir
            else:
                search_path = Path(normalized_dir)
            
            # Glob for all vasprun*.xml files
            matching_files = sorted(search_path.glob("vasprun*.xml"))
            
            if matching_files:
                # Add the normalized paths (relative to base_dir if provided)
                for file_path in matching_files:
                    if base_dir:
                        # Keep only the normalized portion
                        rel_path = str(Path(normalized_dir) / file_path.name)
                    else:
                        rel_path = str(file_path)
                    paths.append(rel_path)
            else:
                # If no files found, add a note or keep the pattern
                # You can choose to either skip or add a warning
                print(f"Warning: No files found matching pattern in {normalized_dir}")
        else:
            # Parse the range and create individual paths
            file_numbers = parse_file_nums(file_nums_spec)
            for num in file_numbers:
                path = str(Path(normalized_dir) / f"vasprun{num}.xml")
                paths.append(path)
    
    return paths


def extract_all_trainval_paths(system: Dict[str, Any], base_dir: Optional[str] = None) -> List[str]:
    """
    Extract all training/validation paths from a system.
    
    Args:
        system: System dictionary containing trainval specifications
        base_dir: Base directory to prepend to normalized paths for globbing
    
    Returns:
        List of normalized file paths
    """
    paths = []
    
    for trainval_entry in system.get("trainval", []):
        paths.extend(extract_trainval_paths(trainval_entry, base_dir))
    
    return paths


def parse_json_file(json_path: str, include_equil: bool = True, base_dir: Optional[str] = None) -> Dict[str, List[str]]:
    """
    Parse a JSON file and extract all vasprun.xml paths organized by system.
    
    Args:
        json_path: Path to the JSON file
        include_equil: Whether to include equilibrium files (False for test sets)
        base_dir: Base directory where the xml files are located (for globbing wildcards)
    
    Returns:
        Dictionary mapping system names to lists of file paths
    """
    with open(json_path, 'r') as f:
        data = json.load(f)
    
    system_paths = {}
    
    for system in data.get("systems", []):
        system_name = system.get("name", "unknown")
        paths = []
        
        # Add equilibrium paths if requested
        if include_equil:
            paths.extend(extract_equil_paths(system))
        
        # Add training/validation paths
        paths.extend(extract_all_trainval_paths(system, base_dir))
        
        system_paths[system_name] = paths
    
    return system_paths


def get_all_paths(json_path: str, include_equil: bool = True, base_dir: Optional[str] = None) -> List[str]:
    """
    Get a flat list of all vasprun.xml paths from a JSON file.
    
    Args:
        json_path: Path to the JSON file
        include_equil: Whether to include equilibrium files
        base_dir: Base directory where the xml files are located (for globbing wildcards)
    
    Returns:
        List of all file paths
    """
    system_paths = parse_json_file(json_path, include_equil, base_dir)
    
    all_paths = []
    for paths in system_paths.values():
        all_paths.extend(paths)
    
    return all_paths


def write_paths_to_file(json_path: str, output_path: str, include_equil: bool = True,
                       organize_by_system: bool = False, base_dir: Optional[str] = None) -> None:
    """
    Parse JSON and write all file paths to a text file.
    
    Args:
        json_path: Path to the JSON file
        output_path: Path to the output text file
        include_equil: Whether to include equilibrium files
        organize_by_system: If True, organize output by system with headers
        base_dir: Base directory where the xml files are located (for globbing wildcards)
                 If None, wildcards will not be expanded and a warning will be printed
    """
    if organize_by_system:
        system_paths = parse_json_file(json_path, include_equil, base_dir)
        
        with open(output_path, 'w') as f:
            for system_name, paths in system_paths.items():
                f.write(f"# System: {system_name}\n")
                for path in paths:
                    f.write(f"{path}\n")
                f.write("\n")
    else:
        all_paths = get_all_paths(json_path, include_equil, base_dir)
        
        with open(output_path, 'w') as f:
            for path in all_paths:
                f.write(f"{path}\n")


def main():
    """Example usage of the module."""
    # Specify the base directory where your xml files are located
    # This should be the directory that contains the xml_collection and xml_set folders
    base_dir = "."  # Change this to your actual data directory
    
    # Example for training/validation set
    trainval_json = "trainval_set_data.json"
    trainval_output = "trainval_paths.txt"
    
    print(f"Processing {trainval_json}...")
    write_paths_to_file(trainval_json, trainval_output, include_equil=True, base_dir=base_dir)
    print(f"Wrote paths to {trainval_output}")
    
    # Example for test set
    test_json = "test_set_data.json"
    test_output = "test_paths.txt"
    
    print(f"\nProcessing {test_json}...")
    write_paths_to_file(test_json, test_output, include_equil=False, base_dir=base_dir)
    print(f"Wrote paths to {test_output}")
    
    # Example: Get paths organized by system
    system_paths = parse_json_file(trainval_json, include_equil=True, base_dir=base_dir)
    print(f"\nFound {len(system_paths)} systems:")
    for system_name, paths in system_paths.items():
        print(f"  {system_name}: {len(paths)} files")


if __name__ == "__main__":
    main()
