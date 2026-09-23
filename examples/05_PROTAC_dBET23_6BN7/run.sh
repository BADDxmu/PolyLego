#!/usr/bin/env bash
# PROTAC ternary complex: CRBN - dBET23 - BRD4(BD1), checked against crystal structure 6BN7.
# Runtime: about 20 minutes (force-field assignment ~8 min, docking ~8 min).
cd "$(dirname "$0")"; source ../common/env.sh
mkdir -p output
TXT=input/dBET23.txt

step "1/5 assemble dBET23 from building blocks (jq1+c8diamine -> part1, part1+thal -> part2)"
PolyLego-BuildPolymer -i input -o output -t $TXT -l part1
PolyLego-BuildPolymer -i input -o output -t $TXT -l part2

step "2/5 assign force-field parameters"
PolyLego-AssignForceField -l input -t $TXT -o output -n part2 -protac --gpu "$GPU"

step "3/5 is the assembled molecule the PROTAC in the crystal (stereocentres included)?"
"$PY" -c "
from PolyLego.protac_dock.dockready import make_dockready_sdf
from polylego_examples import compare_to_crystal_ligand
sdf = make_dockready_sdf('output', name='dBET23')
same, centres = compare_to_crystal_ligand(sdf, 'input/protac_crystal.pdb')
print(sdf, '- same molecule as the crystal ligand:', same, centres)
assert same, 'assembled molecule differs from the crystal ligand'"

step "4/5 ternary-complex poses + MD-ready inputs (pockets from P2Rank)"
# input/dBET23_conformer.sdf is one fixed conformer of dBET23; see examples/README.md for why it is
# fixed and how much the result depends on it (the notebook, section 6).
read -r E3_SITE < <("$PY" -c "from polylego_examples import pocket_center as p; print(*p('input/pockets/e3.pdb_predictions.csv'))")
read -r POI_SITE < <("$PY" -c "from polylego_examples import pocket_center as p; print(*p('input/pockets/poi.pdb_predictions.csv'))")
PolyLego-OptimizeBinding --protacs --md-ready --input_path input/ --output_path output/ \
    --protein1 e3.pdb --protein2 poi.pdb --protac_mol input/dBET23_conformer.sdf \
    --docking_site1 $E3_SITE --docking_site2 $POI_SITE --name dBET23 --topn 5

step "5/5 score the poses against the crystal (DockQ)"
"$PY" -c "
import glob
from polylego_examples import write_reference_complex, dockq, capri_label
ref = write_reference_complex('input/poi.pdb', 'input/e3.pdb', 'output/reference.pdb')
for m in sorted(glob.glob('output/md_inputs/pose_*/complex.pdb')):
    r = dockq(m, ref, {'C'}, {'B'})
    s = r['DockQ'] if r else 0.0
    print(f'{m.split(\"/\")[-2]:10s} DockQ {s:.3f}  ({capri_label(s)})')"
echo; echo "Done: output/md_inputs/pose_00*/complex.pdb, output/reference.pdb"
