import json
from collections import Counter

# Load the JSON file
with open('merged_final_data.json', 'r') as f:
    data = json.load(f)

# Extract all "Three letter code" values
codes = [entry.get("SMILES") for entry in data if "SMILES" in entry]

# Count occurrences
code_counts = Counter(codes)

# Find duplicates
duplicates = {code: count for code, count in code_counts.items() if count > 1}

# Output results
if duplicates:
    print("Duplicate SMILES found:")
    for code, count in duplicates.items():
        print(f"{code}: {count} times")
else:
    print("No duplicate SMILES found.")