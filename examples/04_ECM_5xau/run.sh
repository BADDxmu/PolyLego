#!/usr/bin/env bash
# ECM glycosylation: collagen-binding protein 5XAU (6 chains) + NAG on ASN A3287.
cd "$(dirname "$0")"; source ../common/env.sh
mkdir -p output
step "1/2 protein preparation: ATOM only; truncate disordered A2676-2701 to GLY backbone"
"$PY" -c "
from polylego_examples import strip_hetatm, backbone_only
strip_hetatm('input/5xau.pdb', 'output/5xau_atom.pdb')
backbone_only('output/5xau_atom.pdb', 'output/5xau_protein.pdb', 'A', range(2676, 2702))"
step "2/2 build the conjugate"
PolyLego-BuildConjugate --protein output/5xau_protein.pdb \
    --conjugations input/conjugation_spec.json --output_path output \
    --gmx "$GMX" --gpu "$GPU"
echo; echo "Done: output/conjugate_top/, output/em/"
