#!/usr/bin/env python3
"""
Script to extract and print the mother hierarchy of a geometry from a GDML file.
Usage: python3 gdml_hierarchy.py <gdml_file> <volume_name>
"""

import xml.etree.ElementTree as ET
import sys
from collections import defaultdict

def build_hierarchy(gdml_file):
    """
    Parse GDML file and build a mapping of:
    - volume -> list of children (physvols inside it)
    - volume -> mother (parent)
    """
    tree = ET.parse(gdml_file)
    root = tree.getroot()
    
    # Dictionary to store children of each volume
    children_map = defaultdict(list)
    # Dictionary to store mother of each volume
    mother_map = {}
    # Dictionary to store volume info
    volume_info = {}
    
    # Find structure section
    structure = root.find('.//structure')
    if structure is None:
        print("Error: No structure section found in GDML")
        return None, None, None
    
    # First pass: collect all volumes and physvols
    for volume in structure.findall('volume'):
        vol_name = volume.get('name')
        volume_info[vol_name] = volume
        
        # Find all physical volumes inside this volume
        for physvol in volume.findall('physvol'):
            physvol_name = physvol.get('name')
            volumeref = physvol.find('volumeref')
            if volumeref is not None:
                ref_vol = volumeref.get('ref')
                children_map[vol_name].append((physvol_name, ref_vol))
                mother_map[ref_vol] = vol_name
    
    return children_map, mother_map, volume_info

def print_parents(volume_name, mother_map, indent=0):
    """
    Print the mother hierarchy (all parents up to root) of a volume.
    """
    prefix = "  " * indent
    print(f"{prefix}{volume_name}")
    
    if volume_name in mother_map:
        print_parents(mother_map[volume_name], mother_map, indent + 1)

def print_children(volume_name, children_map, indent=0, max_depth=10):
    """
    Print the children hierarchy (all daughters) of a volume.
    """
    if indent > max_depth:
        return
    
    prefix = "  " * indent
    print(f"{prefix}{volume_name}")
    
    if volume_name in children_map:
        for physvol_name, ref_vol in children_map[volume_name]:
            print(f"{prefix}  └─ {physvol_name} (volume: {ref_vol})")
            print_children(ref_vol, children_map, indent + 2, max_depth)

def search_volume(query, volume_info):
    """
    Search for volumes matching a query (case-insensitive substring match).
    """
    matches = []
    query_lower = query.lower()
    for vol_name in volume_info.keys():
        if query_lower in vol_name.lower():
            matches.append(vol_name)
    return sorted(matches)

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 gdml_hierarchy.py <gdml_file> [volume_name] [--children|--parents|--search]")
        print("")
        print("Options:")
        print("  --children    Print all children (daughters) of the volume (default)")
        print("  --parents     Print all parents (mothers) up to root")
        print("  --search      Search for volumes matching the name pattern")
        sys.exit(1)
    
    gdml_file = sys.argv[1]
    
    try:
        children_map, mother_map, volume_info = build_hierarchy(gdml_file)
    except FileNotFoundError:
        print(f"Error: File not found: {gdml_file}")
        sys.exit(1)
    except ET.ParseError as e:
        print(f"Error parsing GDML file: {e}")
        sys.exit(1)
    
    if children_map is None:
        sys.exit(1)
    
    # Parse arguments
    if len(sys.argv) < 3:
        print("Available volumes (first 20):")
        for i, vol in enumerate(sorted(volume_info.keys())[:20]):
            print(f"  {vol}")
        print(f"  ... and {len(volume_info) - 20} more")
        return
    
    volume_name = sys.argv[2]
    mode = "children"
    
    if len(sys.argv) > 3:
        if sys.argv[3] == "--parents":
            mode = "parents"
        elif sys.argv[3] == "--children":
            mode = "children"
        elif sys.argv[3] == "--search":
            mode = "search"
    
    # Check if volume exists or search for it
    if mode == "search" or volume_name not in volume_info:
        matches = search_volume(volume_name, volume_info)
        if not matches:
            print(f"No volumes found matching: {volume_name}")
            sys.exit(1)
        
        if len(matches) > 1:
            print(f"Found {len(matches)} matching volumes:")
            for match in matches[:20]:
                print(f"  {match}")
            if len(matches) > 20:
                print(f"  ... and {len(matches) - 20} more")
            return
        
        volume_name = matches[0]
    
    print(f"\nVolume: {volume_name}\n")
    
    if mode == "parents":
        print("Mother hierarchy (up to root):")
        print_parents(volume_name, mother_map)
    else:  # children
        print("Children hierarchy (daughters):")
        print_children(volume_name, children_map, max_depth=6)

if __name__ == "__main__":
    main()
