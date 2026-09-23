# PolyLego Tutorial

A walk through what PolyLego does, starting from a three-minute example and ending with
your own system. Every command here is taken from the worked examples in `examples/`,
which are run end to end before each release.

**Contents**

1. [Setup](#1-setup)
2. [Your first assembly](#2-your-first-assembly-3-minutes)
3. [The assembly definition file](#3-the-assembly-definition-file)
4. [Putting a polymer on a protein](#4-putting-a-polymer-on-a-protein)
5. [Covalent conjugates: ADCs and glycosylation](#5-covalent-conjugates-adcs-and-glycosylation)
6. [PROTAC ternary complexes](#6-protac-ternary-complexes)
7. [Continuing to MD](#7-continuing-to-md)
8. [Troubleshooting](#8-troubleshooting)
9. [What to trust](#9-what-to-trust)

---

## 1. Setup

```bash
conda env create -f environment.yml      # or mamba / micromamba
conda activate polylego
pip install PolyLego-1.4.0-cp310-cp310-linux_x86_64.whl
PolyLego-BuildPolymer -h                 # check the install
```

The environment brings its own GROMACS and AmberTools, so nothing else is needed. If you
prefer your own GROMACS build, put `gmx` (or `gmx_mpi`) on `PATH` or set `POLYLEGO_GMX`.

**Threads.** PolyLego runs GROMACS and RDKit underneath. On a many-core machine, pin the
maths libraries or they will fight each other:

```bash
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
```

`mdrun` follows `OMP_NUM_THREADS`; the examples set this for you in `examples/common/env.sh`.

**GPU.** `--gpu auto` (the default) lets GROMACS use a compatible GPU if it finds one, and
the CPU otherwise. `--gpu none` forces CPU; `--gpu 0` pins one device. The same choice is
available as `$POLYLEGO_GPU`.

**Environment variables.**

| Variable | Purpose |
|---|---|
| `POLYLEGO_GMX` | GROMACS executable (default: `gmx` or `gmx_mpi` on `PATH`) |
| `POLYLEGO_GPU` | default for `--gpu`: `auto`, `none`, or a GPU id |
| `POLYLEGO_LOG_LEVEL` | `DEBUG`, `INFO` (default) or `WARNING` |
| `POLYLEGO_FF_CACHE` | cache directory for protein parameters (default `~/.cache/polylego/ff_params`) |

**What a command prints.** A header line, one line per step, and a closing line listing the
files it wrote - all on standard error:

```text
PolyLego 1.4.0 | BuildPolymer
assembling B from input/peg.txt
Finished in 0 s
  output/B.sdf
  output/B_graph.gml
```

Warnings start with `WARNING:`. On failure a command prints one `ERROR:` line and exits
with status 1. Add `-v` for internal diagnostics (and the traceback on failure), `-q` to
keep warnings and errors only. Every command lists its options, grouped, with `-h`; the
full reference is [`CLI_MANUAL.md`](CLI_MANUAL.md).

---

## 2. Your first assembly (3 minutes)

```bash
cd examples/01_polymer_PEG
bash run.sh
```

Two commands run inside:

```bash
PolyLego-BuildPolymer     -i input -o output -t input/peg.txt -l B
PolyLego-AssignForceField -l input -t input/peg.txt -o output -n B --gpu auto
```

* **BuildPolymer** connects the building blocks into the molecular graph named `B`.
  It writes `B.sdf` and `B_graph.gml`. The coordinates in `B.sdf` are a *graph layout*
  (every bond drawn at 1.5 Å) - a picture of the connectivity, not a conformation.
* **AssignForceField** parametrises each building block with GAFF, derives the parameters
  for the junctions between blocks from capped fragments, merges everything into
  `assembly.top`, and relaxes the molecule with a short vacuum EM/NVT. The real 3D
  structure comes out of this step: `assembly.pdb`, `assembly.gro`, `assembly.mol2`.

Check what you got:

```python
import parmed
s = parmed.load_file("output/assembly.top", xyz="output/assembly.gro")
print(len(s.atoms), "atoms,", len(s.bonds), "bonds,",
      sum(d.improper for d in s.dihedrals) + len(s.impropers), "impropers,",
      f"net charge {sum(a.charge for a in s.atoms):+.4f}")
```

For this example: 68 atoms, 67 bonds, 2 impropers, net charge 0.0000.

`assembly.top` + `assembly.gro` are ready for `gmx grompp`.

The notebook `01_polymer_PEG.ipynb` runs the same steps and shows the molecules in 3D.

---

## 3. The assembly definition file

This is the one file you write yourself. It has two tab-separated tables.

```
## [ BUILDING BLOCK ]
name	type	head_smart	head_bonders	head_deletions	end_smart	end_bonder	end_deletor
ethandial	2	[C][O][H]	0	1,2	[C][O][H]	1	2
pei	3	[O]=[C]-[O]-[H]	1	2,3

## [ CONNECTION ]
name	donate_mol	accept_mol	repeat	reaction_times
A	ethandial	ethandial	3	1
B	A	pei	1	1
```

**Rules that bite if you get them wrong:**

* **Columns are read by position, and the header line is skipped.** The names in the header
  are for you, not for the parser - that is why you will see both `head_bonders` and
  `head_bondor` in different examples. What matters is the order:
  `name, type, head_smart, head_bonders, head_deletions, end_smart, end_bonders, end_deletions`.
* **Fields are separated by tabs**, not spaces. A row may stop early (a block with only one
  reactive end writes the first five fields).
* **Every block named anywhere in the file must exist** as `<name>.sdf` in the input
  directory. `AssignForceField` parses the whole table, not only the label you pass with
  `-n`, so a leftover row referring to a missing `.sdf` stops the run.

**`type`** says how many reactive ends the block has:

| `type` | meaning |
|---|---|
| 1 | one reactive end (fill in the `head_*` fields only) |
| 2 | two reactive ends (`head_*` and `end_*`) |
| 3 | one end with several equivalent reactive sites - a branch point or multi-arm core |

**`head_smart` / `end_smart`** is a SMARTS pattern matching the reactive group.
`head_bonders` is the index *within that match* of the atom that forms the new bond, and
`head_deletions` lists the atoms of the match that are removed when it forms (a leaving
group). Indices are 0-based and comma-separated.

Use `ABSOLUTE` instead of a SMARTS pattern to address atoms by their index in the `.sdf`
file directly - this is what the PROTAC example does, because a warhead has exactly one
attachment atom and no pattern is needed:

```
warhead1	1	ABSOLUTE	17	36
```
reads as: atom 17 forms the bond, atom 36 is deleted.

**`[ CONNECTION ]`** builds named assemblies from blocks. `repeat` is how many copies of
`donate_mol` go in; `reaction_times` is how many bonds are formed per connection. An
assembly can be used as a building block of the next row - that is how `A` (three
ethanediols) becomes part of `B` above.

**Stereochemistry.** Stereocentres are taken from the 3D coordinates of each block's
`.sdf` and carried into the assembly, so draw each block as the stereoisomer you mean.
(The PROTAC example checks this: both stereocentres of dBET23 come out as in the crystal.)
A stereocentre created *at* a junction - an atom that only becomes chiral when the new
bond forms - is not defined by the blocks and comes out arbitrary. Assemblies above 500
atoms get their 3D structure from a different generator that does not enforce
stereochemistry; PolyLego warns if a stereocentre came out inverted.

**Block size.** Blocks smaller than ~8 heavy atoms are allowed but warned about: fragment
perception can "see through" a small block from one end to the other, which makes the
junction parameters less specific. The completeness check reports this.

---

## 4. Putting a polymer on a protein

```bash
cd examples/02_polymer_mNET
bash run.sh                       # ~2 minutes
```

After assembling and parametrising the polymer as above:

```bash
PolyLego-OptimizeBinding --polymers --input_path ./ \
    --ligand_pdb assembly.pdb --receptor_pdb mNET.pdb \
    --poccent_file ../input/pockets/mNET.pdb_predictions.csv --name A2
```

`--poccent_file` is a [P2Rank](https://github.com/rdk/p2rank) prediction file
(`prank predict -f protein.pdb`); the examples ship one so you do not need P2Rank
installed. One pose per predicted pocket is written to `output/docking/_<k>/`, as
`A2_<k>_init.pdb` (placed) and `A2_<k>_final_0.3.pdb` (optimised). The pose files contain
the polymer only, in the receptor's coordinate frame - load them together with the
receptor to look at the complex.

---

## 5. Covalent conjugates: ADCs and glycosylation

One command builds a protein with covalently attached modifiers - a drug-linker on an
antibody, a glycan on a glycosylation site:

```bash
cd examples/03_ADC_6p6d
bash run.sh                       # ~2 minutes
```

```bash
PolyLego-BuildConjugate --protein output/6p6d_protein.pdb \
    --conjugations input/conjugation_spec.json --output_path output --gpu auto
```

The spec lists one object per attachment site:

```json
[{"modifier_file": "ligand.sdf", "mol_name": "LIG1",
  "attach_protein_resnum": 297, "attach_protein_atomname": "ND2",
  "attach_protein_occurrence": 1,
  "attach_modifier_idx": 0, "deletor_modifier_idx": 99,
  "deletor_protein_atomname": "HD21"}]
```

* `attach_protein_occurrence` picks *which* residue 297 when several chains have one.
  The ADC example attaches to ASN297 of both heavy chains with `1` and `2`.
* `deletor_*` name the leaving atoms on each side (the placeholder halogen on the modifier,
  the hydrogen on the protein).
* Relative paths inside the JSON are resolved against the JSON file's own directory.

**Protein preparation matters.** `pdb2gmx` only types standard residues:

* strip `HETATM` records (waters, crystallisation additives, glycans you are about to
  rebuild) - `strip_hetatm()` in `examples/common/polylego_examples.py`;
* residues with missing side-chain atoms (disordered in the crystal) stop `pdb2gmx`. The
  glycosylation example (`examples/04_ECM_5xau`) truncates them to a glycine backbone,
  which is a demonstration shortcut - model the missing atoms for production work.

Outputs land in `output/conjugate_top/` (`conjugate.top`, `conjugate.gro`,
`conjugate.pdb` with `LINK` records) and `output/em/` (after vacuum EM). The run prints
the closest heavy-atom distance after EM; check the cross-links themselves:

```python
import parmed, numpy as np
s = parmed.load_file("output/conjugate_top/conjugate.top", xyz="output/em/em/em.gro")
for b in s.bonds:
    r1, r2 = b.atom1.residue.name, b.atom2.residue.name
    if (r1 in {"LIG1", "LIG2"}) != (r2 in {"LIG1", "LIG2"}):
        d = np.linalg.norm(np.array([b.atom1.xx, b.atom1.xy, b.atom1.xz])
                           - np.array([b.atom2.xx, b.atom2.xy, b.atom2.xz]))
        print(f"{r1}-{r2}  {d:.3f} A  (b0 = {b.type.req:.3f} A)")
```

In the ADC example the two ASN-ND2–LIG C1 bonds come out at 1.527 Å and 1.589 Å against
an equilibrium length of 1.522 Å. The exact values move a little from run to run; what
matters is that they sit near `b0` rather than at the 1.82 Å of the pre-positioned
structure written into `conjugate.pdb` before minimisation.

---

## 6. PROTAC ternary complexes

```bash
cd examples/05_PROTAC_dBET23_6BN7
bash run.sh                       # ~20 minutes
```

The PROTAC is assembled from three blocks (warhead + linker + E3 ligand) and parametrised
like any other assembly. Two extra steps are specific to PROTACs:

**Docking-ready molecule.** The assembly step writes bond orders and real geometry into
*different* files (`part2.sdf` has the bond orders but graph-layout coordinates;
`assembly.pdb` has the geometry but no bond orders). Combine them, and record which
building block each atom came from, so the cut points are known rather than guessed:

```python
from PolyLego.protac_dock.dockready import make_dockready_sdf
sdf = make_dockready_sdf("output", name="dBET23")
```

**Ternary docking.** Two modes:

```bash
# candidate poses only
PolyLego-OptimizeBinding --protacs --ensemble  ... --topn 5
# poses + MD-ready input directories + quality control
PolyLego-OptimizeBinding --protacs --md-ready ... --topn 5
```

`--md-ready` writes one directory per pose with `complex.pdb` (both proteins plus the
PROTAC on chain L), `protac.sdf`, `qc.json` and `next_steps.sh`.

**Read `qc.json` before using a pose.** The fields mean:

| field | meaning |
|---|---|
| `min_nb_dist_A` | closest heavy-atom distance between molecules; below 2.0 Å MD will push them apart |
| `n_pairs_below` | how many atom pairs sit under 2.0 / 2.5 / 3.0 Å - one bad contact is different from a whole overlapping patch |
| `warheads_in_pocket` | whether both ends still sit in the pockets they were docked into |
| `em_skipped` / `fmax` / `potential` | energy-minimisation status, if a topology was supplied |
| `pass`, `fail_reasons` | the verdict and, when it fails, why |

**Accuracy.** PolyLego returns a *ranked list of candidates*, not a prediction.
Measured with this release on 68 PROTAC ternary crystal structures from the
PDB (PROTAC conformer generated from SMILES, top 5 poses): 31 could be docked from the
SMILES alone - for the rest the two ends of the PROTAC could not be identified without
building-block information - and in 2 of the 31 the top-ranked pose was CAPRI-acceptable
(DockQ >= 0.23).

The example is a case where it works (top-ranked pose DockQ 0.70, medium) - but only for
the conformer it ships; section 6 of its notebook shows how often other conformers of the
same PROTAC succeed. `examples/06_PROTAC_validation_8G1P` is a case where it does not:
all five candidates are CAPRI-incorrect. Treat the output as candidates to filter with
independent evidence - known binding epitopes, mutagenesis, MD stability.

---

## 7. Continuing to MD

Each MD-ready pose directory carries a `next_steps.sh` with the exact commands:
parametrise the PROTAC (`PolyLego-AssignForceField ... -protac`), then build and simulate
the system (`PolyLego-Analyze --simulation ...`). `PolyLego-Analyze` also does the
analysis side - RMSD, RMSF, radius of gyration, energies, free-energy landscapes and
MM-PBSA. See [`CLI_MANUAL.md`](CLI_MANUAL.md) for its options.

---

## 8. Troubleshooting

| Symptom | Cause and fix |
|---|---|
| `no usable gmx found: set POLYLEGO_GMX to a gmx executable` | Put `gmx`/`gmx_mpi` on `PATH`, or `export POLYLEGO_GMX=/path/to/gmx` |
| `argument -paral: expected one argument` | A value starting with `-` needs the `=` form: `-paral="-ntomp 8 -ntmpi 1"`. Without it, argparse reads `-ntomp` as the `-n` option |
| `Environment variable OMP_NUM_THREADS (8) and the number of threads requested on the command line (4) have different values` | GROMACS refuses the mismatch. PolyLego makes `-ntomp` follow `OMP_NUM_THREADS`; if you pass `-paral` yourself, keep the two consistent |
| `When using GPUs, setting the number of OpenMP threads without specifying the number of ranks...` | Add `-ntmpi 1`, or just use `--gpu` and let PolyLego build the flags |
| `You requested mdrun to use GPU devices with IDs 1, but that includes the following incompatible devices` | That GPU id does not exist or is not usable here - use `--gpu auto` or `--gpu none` |
| `Vacuum EM/NVT of the assembly failed` | The short relaxation inside `AssignForceField` did not run. The message names the log directory; the usual cause is mdrun options that do not fit the machine |
| `<name>Building Block Not Found` (e.g. `pegBuilding Block Not Found`) | A block named in the definition file has no `<name>.sdf` in the input directory - including rows you are not building right now |
| `Residue 'XXX' not found in residue topology database` (pdb2gmx) | Non-standard residues or `HETATM` records reached `pdb2gmx`; strip them first |
| `... does not have the required number of atoms` / missing-atom errors in pdb2gmx | Residues with incomplete side chains; model them, or truncate as the ECM example does |
| `antechamber failed for modifier ...` | Check `output/conjugate_top/modifier_*/antechamber.log`. AmberTools must be on `PATH`; a modifier with unusual valences may need a cleaner input SDF |
| `could not split the PROTAC into anchor fragments` | The PROTAC could not be cut into warhead / linker / E3 ligand. Use `make_dockready_sdf()` so the building-block assignment is recorded instead of guessed |
| `every bond in ... is a single bond` / `every bond in ... is 1.500 A long` | You passed a PDB (no bond orders) or a graph-layout SDF to the docking. Same fix: `make_dockready_sdf()` |
| `ModuleNotFoundError: No module named 'pkg_resources'` | `setuptools` is too new for `openmmforcefields`; the shipped `environment.yml` pins `setuptools=69.5.1` |

---

## 9. What to trust

PolyLego is a modelling tool; different parts of it are supported by different amounts of
evidence, and it is worth knowing which is which.

* **Assembly and force-field generation** are checked on every release: the compiled
  distribution must reproduce the source output field for field, and the topology
  completeness check verifies that every bonded term present in the source fragments
  survives into the assembly.
* **Partial charges are not bit-reproducible.** AM1-BCC charges shift by up to ~0.02 e
  between runs on a loaded machine. Anything you compare across runs should tolerate that.
* **Conjugate junction parameters** are derived per junction from a parametrised fragment;
  the examples check that the cross-link survives energy minimisation at a sensible bond
  length.
* **PROTAC ternary docking is exploratory** - see the numbers in section 6.
* **Nothing here replaces an MD stability check.** A structure that passes the quality
  control still has to survive equilibration.
