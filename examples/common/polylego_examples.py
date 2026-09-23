"""Helpers shared by the PolyLego example notebooks.

Running steps:   check_env(), run(), pocket_center()
Structure prep:  strip_hetatm(), backbone_only()
3D views:        show_mol(), show_grid(), show_complex(), show_gro()

Views use py3Dmol (pip install py3Dmol). They render in Jupyter / JupyterLab;
GitHub's static preview cannot show interactive widgets — open the notebook in
Jupyter or on nbviewer to rotate the molecules.
"""
import os
import shutil
import subprocess
import sys
import time

EXAMPLES = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Same output rules as the PolyLego commands: one handler, known third-party notices
# (e.g. the simtk.openmm deprecation line printed while importing) only with --verbose.
try:
    from PolyLego._log import setup as _setup_logging
    _setup_logging()
except ImportError:          # PolyLego not installed: check_env() reports it
    pass


# ---------------------------------------------------------------------------
# Running steps
# ---------------------------------------------------------------------------

def check_env(need_gmx=True):
    """Print what the examples need and fail early if something is missing."""
    import PolyLego
    print(f"PolyLego   {getattr(PolyLego, '__version__', '?')}  ({os.path.dirname(PolyLego.__file__)})")
    print(f"Python     {sys.version.split()[0]}")
    missing = []
    gmx = os.environ.get("POLYLEGO_GMX") or shutil.which("gmx") or shutil.which("gmx_mpi")
    print(f"GROMACS    {gmx or 'NOT FOUND'}")
    if need_gmx and not gmx:
        missing.append("GROMACS: put gmx (or gmx_mpi) on PATH, or set POLYLEGO_GMX")
    ac = shutil.which("antechamber")
    print(f"AmberTools {ac or 'NOT FOUND'}")
    if not ac:
        missing.append("AmberTools: `conda install -c conda-forge ambertools`")
    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        print(f"{k:22s} {os.environ.get(k, '(unset)')}")
    if missing:
        raise RuntimeError("Missing requirements:\n  - " + "\n  - ".join(missing))
    return gmx


def gmx_binary():
    """The GROMACS executable the examples pass to PolyLego."""
    return os.environ.get("POLYLEGO_GMX") or shutil.which("gmx") or shutil.which("gmx_mpi")


def run(cmd, cwd=None, log=None, tail=15):
    """Run one pipeline step. Shows the last `tail` lines; full output goes to `log`.

    Raises CalledProcessError on a non-zero exit so a failed step never goes unnoticed.
    Returns nothing: in a notebook the CompletedProcess repr would dump the whole
    captured output into the cell a second time.
    """
    cmd = [str(c) for c in cmd]
    print("$ " + " ".join(cmd))
    t0 = time.time()
    p = subprocess.run(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                       universal_newlines=True)
    if log:
        with open(log, "w") as f:
            f.write(p.stdout)
    lines = [l for l in p.stdout.splitlines() if l.strip()]
    for l in lines[-tail:]:
        print("  " + l)
    print(f"  -> exit {p.returncode}, {time.time() - t0:.0f} s")
    if p.returncode != 0:
        raise subprocess.CalledProcessError(p.returncode, cmd, p.stdout)


def pocket_center(predictions_csv, rank=1):
    """(x, y, z) of the rank-th pocket in a P2Rank *_predictions.csv."""
    import csv
    with open(predictions_csv) as f:
        rows = list(csv.DictReader(f, skipinitialspace=True))
    row = {k.strip(): v for k, v in rows[rank - 1].items()}
    return tuple(float(row[k]) for k in ("center_x", "center_y", "center_z"))


# ---------------------------------------------------------------------------
# Structure preparation
# ---------------------------------------------------------------------------

def strip_hetatm(pdb_in, pdb_out):
    """Keep protein ATOM records only (pdb2gmx cannot type glycans/ligands/waters)."""
    with open(pdb_in) as fin, open(pdb_out, "w") as fout:
        for line in fin:
            if not line.startswith("HETATM"):
                fout.write(line)
    return pdb_out


def backbone_only(pdb_in, pdb_out, chain, residues):
    """Truncate the given residues of `chain` to a GLY backbone (N, CA, C, O).

    For disordered crystal residues whose side chains are incomplete: pdb2gmx
    cannot parametrise a residue with missing heavy atoms, but GLY has none.
    """
    residues = {str(r) for r in residues}
    keep = {"N", "CA", "C", "O"}
    with open(pdb_in) as fin, open(pdb_out, "w") as fout:
        for line in fin:
            if line.startswith("ATOM") and line[21] == chain \
                    and line[22:26].strip() in residues:
                if line[12:16].strip() not in keep:
                    continue
                line = line[:17] + "GLY" + line[20:]
            fout.write(line)
    return pdb_out


def gro_to_pdb(top, gro, pdb_out):
    """Write a PDB (with CONECT records) from a GROMACS topology + coordinates."""
    import parmed
    s = parmed.load_file(str(top), xyz=str(gro))
    s.save(str(pdb_out), overwrite=True)
    return pdb_out


# ---------------------------------------------------------------------------
# 3D views
# ---------------------------------------------------------------------------

def _fmt(path):
    ext = os.path.splitext(str(path))[1].lower().lstrip(".")
    return {"sdf": "sdf", "mol": "sdf", "mol2": "mol2", "pdb": "pdb", "xyz": "xyz"}[ext]


def show_mol(path, width=480, height=360, style="stick", label=None):
    """One small molecule / building block / assembled polymer."""
    import py3Dmol
    v = py3Dmol.view(width=width, height=height)
    v.addModel(open(path).read(), _fmt(path))
    v.setStyle({style: {"colorscheme": "Jmol", "radius": 0.18}} if style == "stick"
               else {style: {}})
    if label:
        v.addLabel(label, {"position": {"x": 0, "y": 0, "z": 0}, "fontSize": 12,
                           "backgroundOpacity": 0.6, "inFront": True, "screenOffset": {"x": -200, "y": -150}})
    v.zoomTo()
    return v.show()


def show_grid(paths, titles=None, cols=3, size=260):
    """Several molecules side by side (e.g. the building blocks of an assembly)."""
    import py3Dmol
    rows = (len(paths) + cols - 1) // cols
    v = py3Dmol.view(width=size * cols, height=size * rows, viewergrid=(rows, cols),
                     linked=False)
    for i, p in enumerate(paths):
        r, c = divmod(i, cols)
        v.addModel(open(p).read(), _fmt(p), viewer=(r, c))
        v.setStyle({"stick": {"colorscheme": "Jmol", "radius": 0.18}}, viewer=(r, c))
        if titles:
            v.addLabel(titles[i], {"fontSize": 11, "backgroundOpacity": 0.5,
                                   "useScreen": True, "screenOffset": {"x": 5, "y": 5}},
                       viewer=(r, c))
        v.zoomTo(viewer=(r, c))
    return v.show()


def show_complex(pdb, ligand_resn=(), ligand_chain=None, focus_ligand=True,
                 width=720, height=520, pockets=(), ligands=()):
    """Protein as cartoon (coloured by chain), ligands/modifiers as sticks.

    ligand_resn : residue names drawn as sticks (e.g. ("LIG1", "LIG2"), ("NAG",))
    ligand_chain: or a whole chain drawn as sticks (e.g. "L" for PolyLego PROTAC poses)
    pockets     : (x, y, z) points drawn as translucent spheres
    ligands     : extra structure files (same coordinate frame) drawn as sticks,
                  e.g. a docked polymer pose stored separately from the protein
    """
    import py3Dmol
    v = py3Dmol.view(width=width, height=height)
    v.addModel(open(pdb).read(), "pdb")
    v.setStyle({"cartoon": {"color": "spectrum"}})
    v.setStyle({"chain": list("ABCDEFGHIJ")}, {"cartoon": {"colorscheme": "chain"}})
    sel = []
    if ligand_resn:
        sel.append({"resn": list(ligand_resn)})
    if ligand_chain:
        sel.append({"chain": ligand_chain})
    for s in sel:
        v.setStyle(s, {"stick": {"colorscheme": "greenCarbon", "radius": 0.25}})
        # residues bonded to the ligand, for context
        v.addStyle({"byres": True, "within": {"distance": 4.0, "sel": s}},
                   {"stick": {"radius": 0.12, "colorscheme": "Jmol"}})
    for i, f in enumerate(ligands):
        v.addModel(open(f).read(), _fmt(f))
        v.setStyle({"model": i + 1}, {"stick": {"colorscheme": "magentaCarbon", "radius": 0.25}})
    if ligands and not sel:
        sel.append({"model": 1})
    for (x, y, z) in pockets:
        v.addSphere({"center": {"x": x, "y": y, "z": z}, "radius": 3.0,
                     "color": "orange", "alpha": 0.35})
    if sel and focus_ligand:
        v.zoomTo(sel[0])
        v.zoom(0.35)
    else:
        v.zoomTo()
    return v.show()


def show_gro(top, gro, workdir, **kw):
    """Convenience: GROMACS top+gro -> PDB -> show_complex()."""
    pdb = os.path.join(str(workdir), os.path.splitext(os.path.basename(str(gro)))[0] + "_view.pdb")
    gro_to_pdb(top, gro, pdb)
    return show_complex(pdb, **kw)


# ---------------------------------------------------------------------------
# Validation against a crystal structure
# ---------------------------------------------------------------------------

def smiles_to_sdf(smiles, out_sdf, name="LIG", seed=0xC0FFEE):
    """3D conformer from SMILES (ETKDG + MMFF), heavy atoms only, written as SDF."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    if AllChem.EmbedMolecule(mol, randomSeed=seed) != 0:
        raise RuntimeError("could not embed a 3D conformer")
    AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
    mol = Chem.RemoveHs(mol)
    mol.SetProp("_Name", name)
    w = Chem.SDWriter(str(out_sdf))
    try:
        w.write(mol)
    finally:
        w.close()
    return str(out_sdf)


def write_reference_complex(poi_pdb, e3_pdb, out_pdb):
    """The two chains as they sit in the crystal, written the way PolyLego writes poses.

    Model and reference must come from the same writer so that chain naming and atom
    order match; DockQ itself is superposition-based, so the translation does not matter.
    """
    from PolyLego.protac_dock.geometry import read_pdb, write_complex
    write_complex(out_pdb, read_pdb(poi_pdb), read_pdb(e3_pdb), poi_pdb, e3_pdb)
    return str(out_pdb)


def dockq(model_pdb, reference_pdb, poi_chains, e3_chains, log=None):
    """DockQ of a predicted complex against the crystal (needs the `DockQ` command).

    Returns the best POI-E3 interface as {"DockQ": float, "fnat": ..., "iRMSD": ...},
    or None if DockQ reported no matching interface.
    CAPRI: >= 0.23 acceptable, >= 0.49 medium, >= 0.80 high.
    """
    out = subprocess.run(["DockQ", str(model_pdb), str(reference_pdb),
                          "--allowed_mismatches", "10"],
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                         universal_newlines=True).stdout
    if log:
        open(log, "w").write(out)
    cur, blocks = None, {}
    for line in out.splitlines():
        if line.startswith("Native chains:"):
            cur = tuple(p.strip() for p in line.replace("Native chains:", "").split(","))
            blocks[cur] = {}
        elif cur is not None and "\t" in line and ":" in line:
            k, _, v = line.strip().partition(":")
            try:
                blocks[cur][k.strip()] = float(v.strip())
            except ValueError:
                pass
    hits = [m for (c1, c2), m in blocks.items()
            if ((c1 in poi_chains and c2 in e3_chains)
                or (c1 in e3_chains and c2 in poi_chains)) and "DockQ" in m]
    return max(hits, key=lambda m: m["DockQ"]) if hits else None


def capri_label(score):
    """CAPRI quality class of a DockQ score."""
    return ("high" if score >= 0.80 else "medium" if score >= 0.49
            else "acceptable" if score >= 0.23 else "incorrect")


def compare_to_crystal_ligand(sdf, crystal_pdb):
    """Is the molecule in `sdf` the ligand in `crystal_pdb`, stereocentres included?

    PDB files carry no bond orders, so they are taken from `sdf` and the configuration
    of each stereocentre is read from the crystal coordinates. Returns
    (same: bool, {"assembled": [(atom, CIP)], "crystal": [(atom, CIP)]}).
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem
    mine = Chem.MolFromMolFile(str(sdf))
    Chem.AssignStereochemistryFrom3D(mine)
    from rdkit.rdBase import BlockLogs
    _quiet = BlockLogs()      # symmetric groups give several equivalent matches; harmless here
    xtal = AllChem.AssignBondOrdersFromTemplate(mine, Chem.MolFromPDBFile(str(crystal_pdb)))
    del _quiet
    Chem.AssignStereochemistryFrom3D(xtal)
    centres = {name: Chem.FindMolChiralCenters(m, useLegacyImplementation=False)
               for name, m in (("assembled", mine), ("crystal", xtal))}
    return Chem.MolToSmiles(mine) == Chem.MolToSmiles(xtal), centres
