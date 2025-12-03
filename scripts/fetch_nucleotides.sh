#!/bin/bash
# Fetch the 5 primary nucleotides (bases) for DNA/RNA
# Usage: ./fetch_nucleotides.sh

# List of bases to fetch
BASES=("Adenine" "Guanine" "Cytosine" "Thymine" "Uracil")

echo "Fetching Nucleotides..."

for BASE in "${BASES[@]}"; do
    ./scripts/fetch_molecule.sh "$BASE"
done

echo "------------------------------------------------"
echo "Nucleotide fetch complete."
echo "Note: These are the isolated nitrogenous bases."
echo "To simulate DNA helix formation, you would need to:"
echo "1. Align them via Hydrogen Bonding (A-T = 2 bonds, G-C = 3 bonds)"
echo "2. Connect them via a Phosphate-Deoxyribose backbone."
echo "------------------------------------------------"
