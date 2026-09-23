# PolyLego-BuildConjugate — API & User Manual

**Module version:** 1.1.0
**Packages:** `PolyLego.conjugate_ff_gen`, `PolyLego.conjugate_opt`
**CLI entry point:** `PolyLego-BuildConjugate` → `PolyLego.CLI:build_conjugate`

This manual documents the covalent-conjugation workflow at three levels:
the command line, the JSON spec format, and the Python API.

---

## 1. Concept

`BuildConjugate` builds a GROMACS topology for a protein covalently bonded to
one or more small-molecule modifiers. The pipeline is:

```
protein.pdb ──pdb2gmx──► protein topology ┐
                                           ├─ merge + cross-link ─► conjugate.top/.gro
modifier.sdf ─antechamber► modifier ITP ───┘        │
                                                     ▼
                                            vacuum EM ─► clash check ─► conjugate.pdb (+LINK)
```

Each covalent junction gets **GAFF2 bonded parameters** (bond, angles,
dihedrals, explicit 1-4 pairs) derived from a parametrised junction fragment,
so `grompp` runs without warnings.

---

## 2. Command-line usage

### 2a. Single-site (shorthand, no JSON)

```bash
PolyLego-BuildConjugate \
  --protein         protein.pdb \
  --modifier        drug.sdf \
  --attach_protein  158:SG \
  --attach_modifier 5 \
  --ff_modifier     gaff2 \
  --output_path     ./conjugate_output/
```

### 2b. Multi-site (JSON spec)

```bash
PolyLego-BuildConjugate \
  --protein      protein.pdb \
  --conjugations sites.json \
  --output_path  ./conjugate_output/
```

### 2c. Pre-parametrised modifier (skip antechamber)

```bash
PolyLego-BuildConjugate \
  --protein         protein.pdb \
  --modifier_itp    drug.itp \
  --modifier_gro    drug.gro \
  --attach_protein  158:SG \
  --attach_modifier 5 \
  --output_path     ./conjugate_output/
```

### All options

| Option | Default | Description |
|--------|---------|-------------|
| `--protein` | *(required)* | Protein PDB file |
| `--conjugations` | | Conjugation spec JSON (multi-site) |
| `--output_path` | `./conjugate_output/` | Output directory |
| `--ff_protein` | `amber99sb-ildn` | Protein force field (pdb2gmx) |
| `--water` | `spc` | Water model |
| `--modifier` | | Single modifier SDF/MOL2 (auto-parametrised) |
| `--modifier_itp` | | Pre-parametrised modifier ITP (skips antechamber) |
| `--modifier_gro` | | Modifier GRO (required with `--modifier_itp`) |
| `--attach_protein` | | Attachment site `RESNUM:ATOMNAME`, e.g. `158:SG` |
| `--attach_modifier` | `-1` | Modifier attachment atom index (0-based) |
| `--ff_modifier` | `gaff2` | Modifier force field |
| `--mol_name` | `MOL` | Modifier residue name |
| `--gmx` | `gmx` | GROMACS executable |
| `--parallel_operation` | `-ntomp 4` | GROMACS mdrun flags |
| `--skip_em` | off | Emit topology only, skip vacuum EM |

---

## 3. JSON spec format

`--conjugations` points to a JSON **array**; each object maps to one
`ConjugationSpec`. Multiple objects → multiple sites in a single run.

```json
[
  {
    "modifier_file":            "NAG.sdf",
    "mol_name":                 "NAG",
    "attach_protein_resnum":    297,
    "attach_protein_atomname":  "ND2",
    "attach_modifier_idx":      0,
    "ff_modifier":              "gaff2",
    "attach_protein_occurrence": 0,
    "deletor_modifier_idx":     -1,
    "deletor_protein_atomname": "HD21"
  }
]
```

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `modifier_file` | str | *(req)* | Modifier SDF/MOL2/ITP path |
| `modifier_gro` | str | `""` | GRO (required if `modifier_file` is an ITP) |
| `mol_name` | str | `"MOL"` | Modifier residue name |
| `attach_protein_resnum` | int | *(req)* | Protein residue number (PDB numbering) |
| `attach_protein_atomname` | str | *(req)* | Protein cross-link atom (`ND2`, `SG`, `NZ`, `OG`, …) |
| `attach_modifier_idx` | int | *(req)* | 0-based modifier cross-link atom index |
| `ff_modifier` | str | `"gaff2"` | Modifier force field |
| `attach_protein_occurrence` | int | `0` | Which matching residue when multiple chains share a resnum |
| `deletor_modifier_idx` | int | `-1` | 0-based modifier leaving atom; `-1` = none |
| `deletor_protein_atomname` | str | `None` | Protein leaving atom (e.g. `HD21`); `null` = none |

---

## 4. Python API

### 4a. `PolyLego.conjugate_ff_gen`

```python
from PolyLego.conjugate_ff_gen import (
    ConjugationSpec,
    load_conjugation_specs,
    prepare_protein_topology,
    prepare_modifier_topology,
    build_conjugate_topology,
    find_protein_attachment_atom,
    find_protein_attachment_atom_by_gro,
    find_modifier_attachment_atom,
    apply_single_conjugation,
)
```

| Function | Signature | Returns |
|----------|-----------|---------|
| `ConjugationSpec` | dataclass (fields in §3) | spec object |
| `load_conjugation_specs` | `(json_path)` | `list[ConjugationSpec]` |
| `prepare_protein_topology` | `(protein_pdb, output_dir, ff="amber99sb-ildn", water="spc", gmx="gmx")` | `(top, gro)` |
| `prepare_modifier_topology` | `(modifier_sdf, output_dir, mol_name, ff_type="gaff2")` | `(itp, gro)` |
| `build_conjugate_topology` | `(protein_top, protein_gro, specs, output_dir)` | `(merged_top, merged_gro)` |
| `find_protein_attachment_atom` | `(struct, resnum, atomname, occurrence=0)` | atom index |
| `find_protein_attachment_atom_by_gro` | `(gro_file, struct, resnum, atomname, occurrence=0)` | atom index |
| `find_modifier_attachment_atom` | `(struct, atom_idx)` | atom index |
| `apply_single_conjugation` | `(merged, protein_atom_idx, modifier_atom_global_idx, output_dir)` | cross-link params |

### 4b. `PolyLego.conjugate_opt`

```python
from PolyLego.conjugate_opt import (
    write_vacuum_em_mdp, run_vacuum_em,
    write_vacuum_nvt_mdp, run_vacuum_nvt,
    check_clashes,
)
```

| Function | Signature | Returns |
|----------|-----------|---------|
| `run_vacuum_em` | `(top_file, gro_file, output_dir, gmx="gmx", parallel_options="-ntomp 4", maxwarn=10)` | `em_gro` path |
| `run_vacuum_nvt` | `(top_file, gro_file, output_dir, gmx="gmx", parallel_options="-ntomp 4", n_steps=50000, T=300, maxwarn=10)` | `nvt_gro` path |
| `check_clashes` | `(gro_file, threshold=1.5, top_file=None)` | `(has_clash, min_dist)` |
| `write_vacuum_em_mdp` | `(output_dir)` | mdp path |
| `write_vacuum_nvt_mdp` | `(output_dir, n_steps=50000, T=300, dt=0.001)` | mdp path |

### 4c. End-to-end example

```python
from PolyLego.conjugate_ff_gen import (
    ConjugationSpec, prepare_protein_topology, build_conjugate_topology,
)
from PolyLego.conjugate_opt import run_vacuum_em, check_clashes

# 1. Protein topology
prot_top, prot_gro = prepare_protein_topology(
    "protein.pdb", "./out/protein_top/",
    ff="amber99sb-ildn", water="spc", gmx="gmx",
)

# 2. Define conjugation(s)
specs = [ConjugationSpec(
    modifier_file="NAG.sdf",
    mol_name="NAG",
    attach_protein_resnum=297,
    attach_protein_atomname="ND2",
    attach_modifier_idx=0,
    ff_modifier="gaff2",
    deletor_protein_atomname="HD21",
)]

# 3. Merge + cross-link
merged_top, merged_gro = build_conjugate_topology(
    prot_top, prot_gro, specs, "./out/conjugate_top/",
)

# 4. Vacuum EM
em_gro = run_vacuum_em(merged_top, merged_gro, "./out/em/", gmx="gmx")

# 5. Clash check
has_clash, min_dist = check_clashes(em_gro, top_file=merged_top)
print(f"min heavy-atom distance = {min_dist:.3f} Å  clash={has_clash}")
```

---

## 5. Output files

| File | Description |
|------|-------------|
| `conjugate.top` | GROMACS topology — cross-link bond/angle/dihedral/1-4 pairs |
| `conjugate.gro` | Conjugate structure (post-EM if not `--skip_em`) |
| `conjugate.pdb` | PDB with `LINK` records for ChimeraX / PyMOL |

## 6. Validation

```bash
# Confirm the cross-link bond made it into the topology
grep -A2 "\[ bonds \]" conjugate_output/conjugate.top | head

# grompp should report 0 warnings
gmx grompp -f nvt.mdp -c conjugate_output/conjugate.gro \
           -p conjugate_output/conjugate.top -o nvt.tpr
```

Cross-link chemistry is validated against ADC, ECM, and glycan (Man3GlcNAc2)
systems with 100 ps NVT MD stability checks (tests T7–T10).

---

## 7. Troubleshooting

| Symptom | Cause / fix |
|---------|-------------|
| `grompp` unknown bond type | Junction fragment failed to parametrise — check `--ff_modifier gaff2` and that antechamber (AmberTools) is on `PATH` |
| Clash reported after EM | Modifier pre-positioning collided — try a different `attach_modifier` atom or lower `emtol`; inspect `conjugate.pdb` in ChimeraX |
| Wrong residue bonded (multi-chain) | Set `attach_protein_occurrence` to select the intended chain |
| Extra H on junction | Set `deletor_protein_atomname` / `deletor_modifier_idx` to strip the leaving atom |
