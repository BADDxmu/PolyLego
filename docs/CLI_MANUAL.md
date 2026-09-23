# PolyLego

**Version 1.4.0** · [README](../README.md) · [Tutorial](TUTORIAL.md) · [BuildConjugate API](BuildConjugate_API.md)

PolyLego provides two independent modules for computational chemistry workflows:

| Module | Use case |
|--------|----------|
| **Covalent Conjugation** | Generate GROMACS topologies for covalently modified proteins (N-glycosylation, ADC, etc.) |
| **Polymer–Protein Docking** | PSO-based docking of polymers or PROTAC ternary complexes, with MM-PBSA ΔG |

---

## Output of every command

Each command prints a header line, one line per step, and a closing line with the files
it wrote, all on standard error:

```text
PolyLego 1.4.0 | BuildConjugate
2 conjugation site(s) on 6p6d_protein.pdb
protein topology (pdb2gmx)
modifier parameters and cross-links
vacuum energy minimisation
closest heavy-atom distance after EM: 2.53 A
Finished in 1 min 40 s
  output/conjugate_top/conjugate.top
  output/em/em/em.gro
```

| Option / variable | Effect |
|---|---|
| `-v`, `--verbose` | also print internal diagnostics, with time and module name; show the traceback on failure |
| `-q`, `--quiet` | print warnings and errors only |
| `POLYLEGO_LOG_LEVEL` | `DEBUG`, `INFO` (default) or `WARNING`, for all commands |

Warnings start with `WARNING:`. On failure a command prints one `ERROR:` line and exits
with status 1. Every command lists its options, grouped, with `-h`.

---

## Module 1 — Covalent Conjugation

### Overview

`PolyLego-BuildConjugate` merges a protein topology with one or more small-molecule modifiers at defined covalent attachment sites. It outputs GROMACS-ready `.top` / `.gro` / `.pdb` files with correct cross-link bond, angle, dihedral, and 1-4 non-bonded parameters.

Typical use cases:
- N-linked or O-linked glycosylation (protein–sugar cross-link via ASN ND2 or SER OG)
- Antibody–Drug Conjugates (ADC) — covalent drug attachment to CYS SG or LYS NZ
- Any covalent protein modification requiring explicit GAFF2-parametrised junction chemistry

### Quick Start

```bash
PolyLego-BuildConjugate \
  --protein      protein.pdb \
  --conjugations conjugation_spec.json \
  --output_path  ./conjugate_output/
```

### Conjugation Spec JSON

Define each conjugation site in a JSON file:

```json
[
  {
    "protein_attachment_residue": 297,
    "protein_attachment_atom":   "ND2",
    "protein_deletor_atom":      "HD21",
    "modifier_file":             "NAG.pdb",
    "modifier_attachment_index": 0,
    "modifier_deletor_index":    -1
  }
]
```

| Field | Description |
|-------|-------------|
| `protein_attachment_residue` | Residue number (PDB numbering) |
| `protein_attachment_atom` | Atom on protein forming the cross-link (e.g. `ND2`, `SG`, `NZ`) |
| `protein_deletor_atom` | Hydrogen on protein to remove before bonding (e.g. `HD21`); `null` if none |
| `modifier_file` | Path to modifier PDB (relative to spec file or absolute) |
| `modifier_attachment_index` | 0-based atom index in modifier forming the cross-link bond |
| `modifier_deletor_index` | 0-based index of leaving atom in modifier; `-1` if none |

Multiple entries in the array → multiple conjugation sites processed in one run.

### All Options

| Option | Default | Description |
|--------|---------|-------------|
| `--protein` | *(required)* | Protein PDB file |
| `--conjugations` | | Conjugation spec JSON (multi-site) |
| `--output_path` | `./conjugate_output/` | Output directory |
| `--ff_protein` | `amber99sb-ildn` | Protein force field |
| `--water` | `spc` | Water model |
| `--modifier` | | Single modifier PDB (single-site shorthand, no JSON needed) |
| `--modifier_itp` | | Modifier ITP (pre-parametrised, skips antechamber) |
| `--modifier_gro` | | Modifier GRO |
| `--attach_protein` | | Attachment site as `RESNUM:ATOMNAME`, e.g. `297:ND2` (single-site shorthand) |
| `--attach_modifier` | `-1` | Modifier attachment atom index (0-based) |
| `--ff_modifier` | `gaff2` | Modifier small-molecule force field |
| `--mol_name` | `MOL` | Modifier residue name in ITP |
| `--gmx` | `gmx` | GROMACS executable |
| `--gpu` | `auto` | `auto` (GROMACS picks a compatible GPU, else CPU), `none` (CPU only), or a GPU id such as `0`. Default from `$POLYLEGO_GPU` |
| `--parallel_operation` | from `--gpu` | GROMACS mdrun flags; overrides `--gpu`. Default: `-ntomp $OMP_NUM_THREADS -ntmpi 1` (4 threads if unset) |
| `--skip_em` | | Skip vacuum energy minimisation |

### Output Files

| File | Description |
|------|-------------|
| `conjugate.top` | GROMACS topology — cross-link bond/angle/dihedral/1-4 pairs included |
| `conjugate.gro` | Conjugate structure |
| `conjugate.pdb` | PDB with LINK records for ChimeraX / PyMOL visualisation |

### Validation

After running, confirm the cross-link bond is present in the topology:

```bash
# Check bond entry (1-based atom indices)
grep "RESNUM_ND2_IDX\|MODIFIER_C1_IDX" conjugate_output/conjugate.top

# Validate with GROMACS (expect 0 warnings)
gmx grompp -f nvt.mdp -c conjugate_output/conjugate.gro \
           -p conjugate_output/conjugate.top -o nvt.tpr
```

Visualise in ChimeraX — LINK records are written automatically. If sidechain atoms are hidden by the cartoon representation, run:
```
show :RESNUM atoms
```

---

## Module 2 — Polymer–Protein Docking

### Overview

Four commands form the docking pipeline:

```
PolyLego-BuildPolymer        →  assemble polymer structure
PolyLego-AssignForceField    →  parametrise with GAFF2 / OPLS-AA
PolyLego-OptimizeBinding     →  PSO docking (polymer or PROTAC mode)
PolyLego-Analyze             →  EM + MD + MM-PBSA ΔG
```

### PolyLego-BuildPolymer

Assemble a polymer from building blocks.

```bash
PolyLego-BuildPolymer \
  -i building_block_lib/ \
  -o output_dir/ \
  -t assembly_def.txt \
  -l A2
```

| Option | Description |
|--------|-------------|
| `-i / --input_dir` | Building block structure directory (the `.sdf` files) |
| `-o / --output_dir` | Output directory |
| `-t / --topology_file` | Assembly definition `.txt` — the `[ BUILDING BLOCK ]` + `[ CONNECTION ]` tables, **not** a GROMACS `.top` |
| `-l / --assembly_list` | **Comma**-separated target names from the `[ CONNECTION ]` table (`CLI.py:162` does `.split(",")`), not building block names |

### PolyLego-AssignForceField

Generate force field parameters (GAFF2 or OPLS-AA) for building blocks.

```bash
PolyLego-AssignForceField \
  -l building_block_lib/ \
  -t assembly_def.txt \
  -o output_dir/ \
  -n A2
```

| Option | Description |
|--------|-------------|
| `-l / --bb_lib` | Building block library folder — the same directory passed to `BuildPolymer -i`. Note `-l` means something different here than in `BuildPolymer`. |
| `-t / --top_file` | The same assembly definition `.txt` used by `BuildPolymer` |
| `-o / --output_file` | Output **directory**, despite the name — `CLI.py:205` binds it to `output_dir` and writes ~14 files into it (`assembly.top`, `assembly.txt`, per-fragment `.top`, …) |
| `-n / --bb_name` | Target name, i.e. the same value passed to `BuildPolymer -l` |
| `-protac` | Enable PROTAC mode. **Omit it for polymers** — this flag is the only thing selecting the PROTAC/ternary path. |

Requires **GROMACS `gmx` on `PATH`**: after writing `assembly.top` / `assembly.txt`, `bb_ff_gen.py:1062` shells out to `gmx`. Without it the run ends in `FileNotFoundError: … 'gmx'` *after* force-field generation has already succeeded — so check for the output files before concluding the run failed.

### Choosing CPU or GPU

`PolyLego-AssignForceField`, `PolyLego-Analyze` and `PolyLego-BuildConjugate` all run
GROMACS `mdrun`. Pick the device with `--gpu`:

| `--gpu` | mdrun flags used |
|---|---|
| `auto` (default) | `-ntomp <threads> -ntmpi 1` — GROMACS uses a compatible GPU if it finds one, otherwise the CPU |
| `none` | `... -nb cpu` |
| `0`, `1`, … | `... -nb gpu -gpu_id <id>` |

`<threads>` follows `$OMP_NUM_THREADS` when that is set (GROMACS aborts when the two
disagree). The default can also be set once with `export POLYLEGO_GPU=1`.

To pass mdrun flags yourself, use `-paral` / `--parallel_operation`; it overrides `--gpu`.
**Write it with `=`:**

```bash
PolyLego-AssignForceField ... -paral="-ntomp 8 -ntmpi 1 -nb gpu -gpu_id 0"   # correct
PolyLego-AssignForceField ... -paral "-ntomp 8 -ntmpi 1"                     # also accepted
PolyLego-AssignForceField ... -paral -ntomp 8 -ntmpi 1                       # wrong: -ntomp is read as the -n option
```

The bracket form (`-paral "[-ntomp,8,-ntmpi,1]"`) is still accepted.

### Validation — polymer force-field generation

Verified 2026-07-27 on branch `fix/ff-param-name-collision` @ `577942f`, using the mNET **A2** homopolymer dimer (the "Polymer" row of section 7.4 in `BUGFIX_ff_param_name_collision.md`). Polymer mode = `AssignForceField` **without** `-protac`.

```bash
export PYTHONHASHSEED=0          # makes two runs byte-comparable
PolyLego-BuildPolymer     -i building_blocks/ -o out/ -t A2_only.txt -l A2
PolyLego-AssignForceField -l building_blocks/ -t A2_only.txt -o out/ -n A2 \
                          -paral "[-ntomp,4,-nb,cpu,-ntmpi,1]"

# every atom type must be A_<atomname>, all distinct
awk '/\[ *atoms *\]/{f=1;next} /^\[/{f=0} f&&$1 ~ /^[0-9]+$/{print $2}' out/assembly.top | sort -u
```

Expected for A2:

| Quantity | Value |
|---|---|
| Atom types | 17 distinct, all `A_<atomname>` (`A_C1x`, `A_O1x`, …) |
| `assembly.top` | 31 atoms, 30 bonds, 58 pairs, 102 angles, 207 dihedrals |
| `A.top` (monomer) | 17 atoms, 16 bonds → dimer keeps 29 intra-fragment bonds + 1 cross-link |
| Total charge | −2 × 10⁻⁸ (neutral) |

The two A fragments legitimately share the same 17 type names — that is the repeated-building-block case, and it is correct. Two independent runs in different directories were content-identical; the only differences were the `;   File …` / `;   At date: …` header comments and the absolute `#include ".../assembly.txt"` path, which merely echo the output directory.

Not yet verified: the GROMACS post-processing stage, since `gmx` was absent from the test environment.

> **Known data issue — `data/mNET/A2.txt` cannot be run as shipped.**
> `AssignForceField` parses the *entire* `[ CONNECTION ]` table, not just the label given via `-n`, so it demands every building block named anywhere in the file. `A2.txt` references `H2/H4/H6/H8`, whose `.sdf` files do not exist in the repository — so the recipe in `data/README.md` always fails with `H2Building Block Not Found`. `BuildPolymer` is unaffected, since it only builds the requested label; that is why step 1 succeeds and step 2 does not.
> Workaround: use a trimmed definition file containing only the rows the target needs.

### PolyLego-OptimizeBinding

PSO-based docking optimisation.

```bash
# Polymer–Protein
PolyLego-OptimizeBinding \
  --polymers \
  --protein1    receptor.pdb \
  --head_mol    head.gro \
  --tail_mol    tail.gro \
  --input_path  input/ \
  --output_path output/

# PROTAC ternary complex
PolyLego-OptimizeBinding \
  --protacs \
  --protein1    POI.pdb \
  --protein2    E3.pdb \
  --protac_mol  protac.gro \
  --input_path  input/ \
  --output_path output/
```

| Option | Description |
|--------|-------------|
| `--polymers` | Polymer–protein docking mode |
| `--protacs` | PROTAC ternary complex mode |
| `--protein1` | Primary receptor PDB |
| `--protein2` | Secondary receptor PDB (PROTAC mode) |
| `--head_mol` | Polymer head GRO |
| `--tail_mol` | Polymer tail GRO |
| `--protac_mol` | PROTAC molecule GRO |
| `--input_path` | Input directory |
| `--output_path` | Output directory |
| `--grid_size` | Docking grid size (default `2.0`) |

### PolyLego-Analyze

Run GROMACS EM + NVT/NPT/MD and MM-PBSA analysis.

```bash
PolyLego-Analyze \
  -simulation -simple \
  -pdb    complex.pdb \
  -lig_gro ligand.gro \
  -lig_itp ligand.itp \
  -em_mdp  em.mdp -nvt_mdp nvt.mdp -npt_mdp npt.mdp -md_mdp md.mdp \
  -mm
```

| Option | Description |
|--------|-------------|
| `-simulation` | Run full simulation pipeline |
| `-simple` | Single-ligand mode |
| `-pair` | Pairwise binding mode |
| `-trio` | Ternary complex mode |
| `-vac` | Vacuum simulation |
| `-pdb` | Input PDB |
| `-lig_gro` | Ligand GRO |
| `-lig_itp` | Ligand ITP |
| `-em_mdp / -nvt_mdp / -npt_mdp / -md_mdp` | MDP parameter files |
| `-mm` | Run MM-PBSA analysis |
| `-rmsd / -rmsf / -energy` | Additional trajectory analyses |
| `-ff` | Force field (default `oplsaa`) |
| `--gpu` | GPU for mdrun: `auto` (default), `none`, or a GPU id; also `$POLYLEGO_GPU` |
| `-paral` | GROMACS mdrun flags; overrides `--gpu`. Write it as `-paral="-ntomp 8 -ntmpi 1"` (see note below). Default: `-ntomp $OMP_NUM_THREADS -ntmpi 2` (8 threads if unset) |

### Full Pipeline Example

```bash
# 1. Assemble polymer   (-t is the assembly definition .txt; -l is a name from its [ CONNECTION ] table)
PolyLego-BuildPolymer -i bb_lib/ -o polymer/ -t peg_def.txt -l A

# 2. Parametrise        (no -protac == polymer path; -o is a directory)
PolyLego-AssignForceField -l bb_lib/ -t peg_def.txt -o polymer/ -n A

# 3. Dock
PolyLego-OptimizeBinding \
  --polymers --protein1 receptor.pdb \
  --head_mol polymer/head.gro --tail_mol polymer/tail.gro \
  --input_path input/ --output_path docking/

# 4. Simulate and score
PolyLego-Analyze -simulation -simple \
  -pdb docking/complex.pdb \
  -lig_gro polymer/polymer.gro -lig_itp polymer/ligand.itp \
  -em_mdp mdp/em.mdp -nvt_mdp mdp/nvt.mdp -npt_mdp mdp/npt.mdp -md_mdp mdp/md.mdp \
  -mm
```

---

## Installation

### From source (development, editable)

```bash
conda activate polylego
pip install -e .
```

### From a built distribution

```bash
conda activate polylego

# Pure-Python wheel (source visible)
pip install dist/PolyLego-1.1.0-py3-none-any.whl

# Cython-compiled wheel (source-protected, Linux / CPython 3.10)
pip install dist/polylego_compiled-1.1.0-cp310-cp310-linux_x86_64.whl
```

### Building the distributions yourself

```bash
conda activate polylego

# Pure-Python sdist + wheel  →  dist/PolyLego-1.1.0-*
python setup.py sdist bdist_wheel

# Cython-compiled wheel  →  dist/polylego_compiled-1.1.0-*
python build_scripts/setup_cython.py bdist_wheel
```

**Dependencies:** GROMACS, p2rank, RDKit, BioPython, parmed, antechamber (AmberTools)

**Repository:** https://gitee.com/JoanDu/PolyLego
