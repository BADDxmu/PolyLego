#!/usr/bin/env bash
# Polymer: PEG-like chain (3 x ethanediol + a PEI cap), 68 atoms. ~3 min.
cd "$(dirname "$0")"; source ../common/env.sh
mkdir -p output
step "1/2 assemble the polymer from building blocks"
PolyLego-BuildPolymer -i input -o output -t input/peg.txt -l B
step "2/2 assign force-field parameters (GAFF building blocks + junction parameters)"
PolyLego-AssignForceField -l input -t input/peg.txt -o output -n B --gpu "$GPU"
echo; echo "Done: output/assembly.pdb, output/assembly.top"
