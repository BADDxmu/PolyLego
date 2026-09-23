#!/usr/bin/env bash
# Polymer-protein: A2 dimer docked onto the mNET protein.
cd "$(dirname "$0")"; source ../common/env.sh
mkdir -p output
step "1/3 assemble the A2 polymer"
PolyLego-BuildPolymer -i input -o output -t input/A2.txt -l A2
step "2/3 assign force-field parameters"
PolyLego-AssignForceField -l input -t input/A2.txt -o output -n A2 --gpu "$GPU"
step "3/3 optimise binding of A2 on mNET (pocket from P2Rank)"
cp input/mNET.pdb output/
( cd output && PolyLego-OptimizeBinding --polymers --input_path ./ \
      --ligand_pdb assembly.pdb --receptor_pdb mNET.pdb \
      --poccent_file ../input/pockets/mNET.pdb_predictions.csv --name A2 )
echo; echo "Done: see output/"
