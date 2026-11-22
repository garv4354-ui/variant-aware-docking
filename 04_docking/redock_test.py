import subprocess, os, numpy as np
from scipy.optimize import linear_sum_assignment

CENTER = (22.01, 0.25, 52.79)
os.makedirs("redock_test", exist_ok=True)

def read_models(path):
    models, cur, scores = [], {}, []
    for l in open(path):
        if l.startswith("MODEL"): cur = {}
        elif l.startswith("REMARK VINA RESULT"): scores.append(float(l.split()[3]))
        elif l.startswith("ENDMDL"): models.append({k: np.array(v) for k, v in cur.items()})
        elif l.startswith(("ATOM", "HETATM")):
            el = l[12:16].strip()[0].upper()
            if el == "H": continue
            cur.setdefault(el, []).append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return models, scores

def approx_rmsd(a, b):
    tot = n = 0
    for el in a:
        D = np.linalg.norm(a[el][:, None] - b[el][None], axis=2)
        r, c = linear_sum_assignment(D)
        tot += (D[r, c] ** 2).sum(); n += len(r)
    return (tot / n) ** 0.5

crys, _ = read_models("../03_structures/prepared/erlotinib_crystal.pdb")
if not crys:   # crystal file has no MODEL records: read it as one model
    d = {}
    for l in open("../03_structures/prepared/erlotinib_crystal.pdb"):
        if l.startswith("HETATM"):
            el = l[12:16].strip()[0].upper()
            if el == "H": continue
            d.setdefault(el, []).append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    crys = [{k: np.array(v) for k, v in d.items()}]
crys = crys[0]

for box, exh in [(22, 16), (14, 32)]:
    out = f"redock_test/erlotinib_box{box}_exh{exh}.pdbqt"
    subprocess.run(["vina", "--receptor", "receptors/EGFR_WT.pdbqt",
        "--ligand", "ligands/erlotinib.pdbqt",
        "--center_x", str(CENTER[0]), "--center_y", str(CENTER[1]), "--center_z", str(CENTER[2]),
        "--size_x", str(box), "--size_y", str(box), "--size_z", str(box),
        "--exhaustiveness", str(exh), "--num_modes", "9", "--energy_range", "4",
        "--seed", "1", "--cpu", "2", "--out", out], check=True, capture_output=True)
    models, scores = read_models(out)
    print(f"\nbox {box} A, exhaustiveness {exh}")
    for i, (m, s) in enumerate(zip(models, scores), 1):
        print(f"  pose {i}: score {s:6.2f}   approx RMSD {approx_rmsd(crys, m):5.2f} A")
