#!/usr/bin/env bash
# ADC: IgG1 Fc (PDB 6P6D, 2 chains) with a drug-linker on ASN297 of each chain.
cd "$(dirname "$0")"; source ../common/env.sh
mkdir -p output
step "1/2 protein preparation: keep ATOM records only"
"$PY" -c "from polylego_examples import strip_hetatm; strip_hetatm('input/6p6d.pdb', 'output/6p6d_protein.pdb')"
step "2/2 build the conjugate: pdb2gmx + modifier GAFF2 + cross-links + vacuum EM"
PolyLego-BuildConjugate --protein output/6p6d_protein.pdb \
    --conjugations input/conjugation_spec.json --output_path output \
    --gmx "$GMX" --gpu "$GPU"
echo; echo "Done: output/conjugate_top/, output/em/"
