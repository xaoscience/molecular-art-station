#!/bin/bash
set -e

echo "--- Installing ChemGen Dependencies ---"
echo "Updating apt..."
sudo apt-get update

echo "Installing C++ Build Tools and SDL2..."
sudo apt-get install -y build-essential libsdl2-dev libsdl2-ttf-dev curl gawk

echo "Done! You can now run ./build.sh"
