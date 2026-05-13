import json

def count_non_nitrogen_entries(filename):
    try:
        with open(filename, 'r') as f:
            data = json.load(f)
    except FileNotFoundError:
        print(f"Error: {filename} not found.")
        return

    no_nitrogen_count = 0
    no_nitrogen_list = []

    for entry in data:
        smiles = entry.get("SMILES", "")
        
        # Check if 'N' or 'n' is NOT in the SMILES string
        if 'N' not in smiles and 'n' not in smiles:
            no_nitrogen_count += 1
            no_nitrogen_list.append(entry.get("Actual PDB code", "Unknown PDB"))

    print(f"--- Analysis of {filename} ---")
    print(f"Number of entries with no Nitrogen ('N') in SMILES: {no_nitrogen_count}")
    
    if no_nitrogen_list:
        print(f"PDB Codes found: {', '.join(no_nitrogen_list)}")

if __name__ == "__main__":
    count_non_nitrogen_entries('merged_final_data.json')