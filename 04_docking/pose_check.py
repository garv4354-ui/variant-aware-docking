import sys, numpy as np
from scipy.optimize import linear_sum_assignment

def read(path):
    d = {}
    for l in open(path):
        if l.startswith("ENDMDL"): break
        if l.startswith(("ATOM", "HETATM")):
            el = l[12:16].strip()[0].upper()
            if el == "H": continue
            d.setdefault(el, []).append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return {k: np.array(v) for k, v in d.items()}

crys = read("../03_structures/prepared/erlotinib_crystal.pdb")
dock = read(sys.argv[1] if len(sys.argv) > 1 else "results/EGFR_WT__erlotinib.pdbqt")
print("atoms crystal:", {k: len(v) for k, v in crys.items()})
print("atoms docked: ", {k: len(v) for k, v in dock.items()})
tot = n = 0
for el in crys:
    D = np.linalg.norm(crys[el][:, None] - dock[el][None], axis=2)
    r, c = linear_sum_assignment(D)
    tot += (D[r, c] ** 2).sum(); n += len(r)
print("approx RMSD (A):", round(float(np.sqrt(tot / n)), 2))
allc = np.vstack(list(crys.values())); alld = np.vstack(list(dock.values()))
print("centroid distance (A):", round(float(np.linalg.norm(allc.mean(0) - alld.mean(0))), 2))
