"""How much does the docking result of this example depend on the input conformer?

Reruns the comparison in section 6 of the notebook:
  * 9 RDKit conformers of dBET23 (seeds 1-8 and 0xC0FFEE), docked straight from memory;
  * the shipped conformer with every atom moved by ~1e-4 A (8 random draws);
  * all poses of the 9 conformers pooled and ranked with the same score.

    python robustness.py            # about 15 minutes with 17 cores free

Run it after run.sh (it needs output/reference.pdb). Conformers are docked from memory,
not written to SDF first: rounding the coordinates to the 4 decimals of an SDF file is
itself enough to change the result, which is the point of this script.
"""
import json
import os
import sys
from multiprocessing import Pool
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS",
           "MMFF_THREADS", "POLYLEGO_MMFF_THREADS"):
    os.environ[_v] = "1"

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "common"))
OUT = HERE / "output" / "robustness"


def _job(spec):
    import numpy as np
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Geometry import Point3D
    from polylego_examples import pocket_center, dockq
    from PolyLego.protac_dock import dock_ternary
    from PolyLego.protac_dock.geometry import read_pdb, write_complex

    kind, seed = spec
    if kind == "conformer":
        smiles = (HERE / "input/dBET23.smi").read_text().split()[0]
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        if AllChem.EmbedMolecule(mol, randomSeed=seed) != 0:
            return {"kind": kind, "seed": seed, "error": "embedding failed"}
        AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
        mol = Chem.RemoveHs(mol)
    else:                                              # "perturbed"
        mol = Chem.MolFromMolFile(str(HERE / "input/dBET23_conformer.sdf"))
        rng, conf = np.random.default_rng(seed), mol.GetConformer()
        for i in range(mol.GetNumAtoms()):
            p, d = conf.GetAtomPosition(i), rng.normal(scale=1e-4, size=3)
            conf.SetAtomPosition(i, Point3D(p.x + d[0], p.y + d[1], p.z + d[2]))
    e3_site = list(pocket_center(str(HERE / "input/pockets/e3.pdb_predictions.csv")))
    poi_site = list(pocket_center(str(HERE / "input/pockets/poi.pdb_predictions.csv")))
    poi_pdb, e3_pdb = HERE / "input/poi.pdb", HERE / "input/e3.pdb"
    poses = dock_ternary(poi_pdb, e3_pdb, mol, poi_site, e3_site, topn=5)
    e3 = read_pdb(e3_pdb)
    d = OUT / f"{kind}_{seed}"
    d.mkdir(parents=True, exist_ok=True)
    rows = []
    for r in poses:
        f = d / f"pose_{r['rank']}.pdb"
        write_complex(f, r["poi_xyz"], e3, poi_pdb, e3_pdb)
        q = dockq(str(f), str(HERE / "output/reference.pdb"), {"C"}, {"B"})
        rows.append({"rank": r["rank"], "dockq": round(q["DockQ"], 3) if q else 0.0,
                     "lj": r["lj"], "n_contact": r["n_contact"], "n_clash": r["n_clash"]})
    return {"kind": kind, "seed": seed, "poses": rows}


def main():
    if not (HERE / "output/reference.pdb").exists():
        sys.exit("run run.sh first (output/reference.pdb is missing)")
    import numpy as np
    specs = [("conformer", s) for s in (1, 2, 3, 4, 5, 6, 7, 8, 0xC0FFEE)]
    specs += [("perturbed", s) for s in range(1, 9)]
    with Pool(min(len(specs), os.cpu_count() or 1)) as pool:
        res = pool.map(_job, specs)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "results.json").write_text(json.dumps(res, indent=1))

    def summary(kind):
        rs = [r for r in res if r["kind"] == kind and "poses" in r]
        ok = sum(r["poses"][0]["dockq"] >= 0.23 for r in rs)
        print(f"{kind:10s} top-ranked pose acceptable in {ok} of {len(rs)}")
        for r in rs:
            print(f"    seed {r['seed']:>9}  top-1 DockQ {r['poses'][0]['dockq']:.3f}"
                  f"   all {[p['dockq'] for p in r['poses']]}")
    summary("conformer")
    summary("perturbed")

    # pooled ranking: the same key dock_ternary uses inside one run
    P = [p for r in res if r["kind"] == "conformer" and "poses" in r for p in r["poses"]]
    ncs = -np.array([p["n_contact"] for p in P], float)
    lj = np.array([p["lj"] for p in P])
    key = (ncs - ncs.mean()) / (ncs.std() + 1e-9) + (lj - lj.mean()) / (lj.std() + 1e-9)
    from PolyLego.protac_dock.pipeline import DEFAULTS
    clashes = np.array([p["n_clash"] for p in P])
    key = np.where(clashes <= DEFAULTS["max_clash"], key, key.max() + 1.0 + clashes)
    order = list(np.argsort(key))
    best = max(range(len(P)), key=lambda i: P[i]["dockq"])
    print(f"pooled     best pose (DockQ {P[best]['dockq']:.3f}) ranks "
          f"{order.index(best) + 1} of {len(P)}")


if __name__ == "__main__":
    main()
