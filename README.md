# ChemGen: Molecular Art Station

A standalone, low-level generative art engine that uses chemical rules and NCBI data to create dynamic, scientifically-inspired visual effects.

## Concept
ChemGen is a C++ based particle engine where particles behave like atoms. They respect valency rules (e.g., Carbon forms 4 bonds), form dynamic covalent bonds based on proximity, and can be influenced by user interaction (mouse cursor).

It bridges "Real Science" and "Generative Art" by fetching real molecular structures from the NCBI PubChem database using Bash scripts, converting them into templates, and then using those templates to spawn "chemical debris" that reassembles procedurally.

First demo:
https://youtu.be/3W1T1ZDGFEk

## Tech Stack
- **Language:** C++ (No Python, no heavy frameworks)
- **Graphics:** SDL2 (Simple DirectMedia Layer) for low-level windowing and rendering.
- **Data:** Bash + Awk scripts to fetch and parse SDF files from NCBI PubChem.
- **Build:** G++ / Make.

## Prerequisites
- Linux (Ubuntu/Debian)
- `apt` package manager
- `g++`, `make`
- `libsdl2-dev`, `libsdl2-ttf-dev`

## Quick Start

### 1. Install Dependencies
```bash
./install_deps.sh
```

### 2. Fetch Data (Optional)
Downloads real molecular structures (e.g., Caffeine, Benzene) from NCBI to use as spawn templates.
```bash
./scripts/fetch_molecule.sh caffeine
./scripts/fetch_molecule.sh benzene
```

### 3. Build
```bash
./build.sh
```

### 4. Run
```bash
./chemgen
```

## Controls
- **Mouse Move:** Emits atoms/radicals.
- **Click:** Spawns a complex molecule from the fetched templates.
- **Space:** Clear screen.
- **ESC:** Quit.
