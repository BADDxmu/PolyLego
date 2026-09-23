#!/usr/bin/env bash
# Validation: predict the SMARCA2 - PROTAC - VHL ternary complex of PDB 8G1P from
# the two proteins and the PROTAC SMILES, then score the poses against the crystal
# with DockQ. Runtime: about 10 minutes.
cd "$(dirname "$0")"; source ../common/env.sh
mkdir -p output

step "1/3 PROTAC 3D conformer from SMILES (no crystal coordinates)"
"$PY" -c "
from polylego_examples import smiles_to_sdf
smiles = open('input/protac.smi').read().split()[0]
print(smiles_to_sdf(smiles, 'output/protac.sdf', name='FWZ'))"

step "2/3 ternary docking (binding sites from P2Rank)"
read -r E3_SITE < <("$PY" -c "from polylego_examples import pocket_center as p; print(*p('input/pockets/e3.pdb_predictions.csv'))")
read -r POI_SITE < <("$PY" -c "from polylego_examples import pocket_center as p; print(*p('input/pockets/poi.pdb_predictions.csv'))")
PolyLego-OptimizeBinding --protacs --ensemble \
    --input_path input/ --output_path output/ \
    --protein1 e3.pdb --protein2 poi.pdb --protac_mol output/protac.sdf \
    --docking_site1 $E3_SITE --docking_site2 $POI_SITE \
    --name 8G1P --topn 5

step "3/3 score the poses against the crystal (DockQ)"
"$PY" -c "
import glob
from polylego_examples import write_reference_complex, dockq, capri_label
ref = write_reference_complex('input/poi.pdb', 'input/e3.pdb', 'output/reference.pdb')
for m in sorted(glob.glob('output/docking_ensemble/*.pdb')):
    r = dockq(m, ref, {'G'}, {'A', 'B', 'C', 'E'})
    s = r['DockQ'] if r else 0.0
    print(f'{m.split(\"/\")[-1]:16s} DockQ {s:.3f}  ({capri_label(s)})')"
echo; echo "Done: output/docking_ensemble/, output/reference.pdb"
