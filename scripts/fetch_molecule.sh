#!/bin/bash
# Fetch a molecule from NCBI PubChem in SDF format and parse it into a simple .chem format for our C++ engine
# Usage: ./fetch_molecule.sh <name>

NAME=${1:-benzene}
OUT_DIR="data"
mkdir -p "$OUT_DIR"
OUT_FILE="$OUT_DIR/$NAME.chem"

echo "Fetching structure for '$NAME' from NCBI PubChem..."

# Download SDF (Structure Data File)
# SDF format:
# Header
# Counts line (aaabbb...) where aaa is num atoms, bbb is num bonds
# Atom block (x y z symbol...)
# Bond block (1 2 type...)

RAW_SDF=$(curl -s "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/$NAME/SDF")

if [[ -z "$RAW_SDF" ]]; then
    echo "Error: Could not fetch data for $NAME"
    exit 1
fi

echo "Parsing SDF..."

# Use awk to extract atoms and bonds
# We create a simple format:
# ATOM <id> <x> <y> <element>
# BOND <id1> <id2> <type>

echo "$RAW_SDF" | awk '
    BEGIN { section = "header"; atom_count = 0; bond_count = 0; atoms_read = 0; bonds_read = 0; }
    
    # The counts line is usually the 4th line in SDF, but it is tricky to identify strictly by line number if header varies.
    # Standard SDF: Line 4 has "aaabbb..." format.
    NR == 4 {
        atom_count = substr($0, 1, 3) + 0
        bond_count = substr($0, 4, 3) + 0
        section = "atoms"
        print "# Molecule: " name
        print "COUNTS " atom_count " " bond_count
        next
    }

    section == "atoms" && atoms_read < atom_count {
        # SDF Atom Line: x(10) y(10) z(10) symbol(3) ...
        x = substr($0, 1, 10) + 0
        y = substr($0, 11, 10) + 0
        symbol = substr($0, 32, 3)
        gsub(/ /, "", symbol) # Trim spaces
        
        print "ATOM " (atoms_read + 1) " " x " " y " " symbol
        atoms_read++
        
        if (atoms_read == atom_count) {
            section = "bonds"
        }
        next
    }

    section == "bonds" && bonds_read < bond_count {
        # SDF Bond Line: 1(3) 2(3) type(3) ...
        a1 = substr($0, 1, 3) + 0
        a2 = substr($0, 4, 3) + 0
        type = substr($0, 7, 3) + 0
        
        print "BOND " a1 " " a2 " " type
        bonds_read++
        next
    }
' name="$NAME" > "$OUT_FILE"

echo "Saved template to $OUT_FILE"
