# PolyLego examples

| Example | Demonstrates | Main commands | Runtime* |
|---|---|---|---|
| `01_polymer_PEG` | building-block assembly and force-field assignment | BuildPolymer, AssignForceField | ~3 min |
| `02_polymer_mNET` | polymer binding on a protein surface | + OptimizeBinding `--polymers` | ~2 min |
| `03_ADC_6p6d` | two covalent drug-linker conjugations on an antibody Fc | BuildConjugate | ~2 min |
| `04_ECM_5xau` | glycan (NAG) attachment on a multi-chain protein | BuildConjugate | ~2 min |
| `05_PROTAC_dBET23_6BN7` | PROTAC assembly, ternary-complex poses and MD-ready inputs, scored against crystal 6BN7 | BuildPolymer, AssignForceField, OptimizeBinding `--protacs --md-ready` | ~20 min |
| `06_PROTAC_validation_8G1P` | the prediction from SMILES only, scored against crystal 8G1P | OptimizeBinding `--protacs --ensemble` | ~10 min |

\* 8 threads, measured on a Linux workstation.

## Running

Every example can be run two ways, from its own directory:

```bash
cd examples/01_polymer_PEG
bash run.sh                       # command line; results in output/
jupyter lab 01_polymer_PEG.ipynb  # notebook: same steps plus 3D views
```

The notebooks use [py3Dmol](https://pypi.org/project/py3Dmol/) for the 3D views. The
committed notebooks contain the text output of a reference run; the 3D views appear when
you run them in Jupyter.

GPU use for the GROMACS runs is chosen with `GPU`: `auto` (default - GROMACS uses a
compatible GPU if it finds one, otherwise the CPU), `none` (CPU only), or a GPU id:
`GPU=0 bash run.sh`. In the notebooks, set `GPU` in the setup cell. The same choice is the
`--gpu` option of `PolyLego-AssignForceField`, `PolyLego-BuildConjugate` and
`PolyLego-Analyze` (or the `POLYLEGO_GPU` environment variable).

`common/` holds the helpers the notebooks import (`polylego_examples.py`) and the shared
shell settings for `run.sh` (`env.sh`).

## Notes on the inputs

* Protein structures are from the RCSB PDB (6P6D, 5XAU, 6BN7, 8G1P); crystal waters and
  hetero groups are removed where the example says so.
* Pocket centres (`input/pockets/*_predictions.csv`) are P2Rank predictions shipped so the
  examples do not need P2Rank. Re-create them with `prank predict -f <protein>.pdb`.
* `04_ECM_5xau` truncates the disordered residues A2676-2701 to a glycine backbone so that
  `pdb2gmx` can type them. This is a shortcut for the demonstration; model the missing atoms
  for production work.
* `05_PROTAC_dBET23_6BN7` docks one fixed conformer of dBET23, `input/dBET23_conformer.sdf`,
  and its top-ranked pose is medium quality against the crystal (DockQ 0.70). That conformer
  was chosen because it works: with 9 other conformers of the same molecule the top pose
  was acceptable once, and moving every atom of the shipped one by ~0.0001 A keeps the
  result in 3 of 8 tries. Section 6 of the notebook has the numbers and `robustness.py`
  reruns them. Identical inputs give identical poses, across runs and thread counts.
* `06_PROTAC_validation_8G1P` starts from the PROTAC SMILES only; on that case all five
  poses are CAPRI-incorrect.
* Measured with this release on 68 PROTAC ternary crystal structures from the
  PDB (PROTAC conformer generated from SMILES, top 5 poses): 31 could be docked from the
  SMILES alone - for the rest the two ends of the PROTAC could not be identified without
  building-block information - and in 2 of the 31 the top-ranked pose was CAPRI-acceptable
  (DockQ >= 0.23). Treat the poses as candidates to be checked
  against other evidence.
