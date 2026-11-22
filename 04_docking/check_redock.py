import numpy as np
from scipy.optimize import linear_sum_assignment

def models(path):
    out, cur = [], {}
    for l in open(path):
        if l.startswith("MODEL"): cur = {}
        elif l.startswith("ENDMDL"): out.append({k: np.array(v) for k, v in cur.items()})
        elif l.startswith(("ATOM", "HETATM")):
            el = l[12:16].strip()[0].upper()
            if el != "H":
                cur.setdefault(el, []).append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return out

def rmsd(a, b):
    t = n = 0
    for el in a:
        D = np.linalg.norm(a[el][:, None] - b[el][None], axis=2)
        r, c = linear_sum_assignment(D)
        t += (D[r, c] ** 2).sum(); n += len(r)
    return (t / n) ** 0.5

d = {}
for l in open("../03_structures/prepared/erlotinib_crystal.pdb"):
    if l.startswith("HETATM"):
        el = l[12:16].strip()[0].upper()
        if el != "H":
            d.setdefault(el, []).append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
crys = {k: np.array(v) for k, v in d.items()}

for f in ["redock_test/erlotinib_box22_exh16.pdbqt", "redock_test/erlotinib_box14_exh32.pdbqt"]:
    r = [round(rmsd(crys, m), 2) for m in models(f)]
    print(f.split("/")[1], "poses:", len(r), "RMSDs:", r, "best:", min(r))
