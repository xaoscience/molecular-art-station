#!/bin/bash
set -e

# Compile the C++ engine
echo "Compiling ChemGen..."

g++ -o chemgen src/main.cpp -lSDL2 -lSDL2_ttf -I./src

echo "Build complete. Run ./chemgen to start."
