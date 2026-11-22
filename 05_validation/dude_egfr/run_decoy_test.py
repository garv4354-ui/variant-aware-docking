import os, random, csv, subprocess
from rdkit import Chem
from rdkit.Chem import AllChem

random.seed(7)
N_ACT, N_DEC = 20, 80
REC = "../../04_docking/receptors/EGFR_WT.pdbqt"
CENTER = (22.01, 0.25, 52.79)
os.makedirs("lig", exist_ok=True); os.makedirs("out", exist_ok=True)

def read_ism(path):
    return [l.split()[0] for l in open(path) if l.strip()]

acts = random.sample(read_ism("actives_final.ism"), N_ACT)
decs = random.sample(read_ism("decoys_final.ism"), N_DEC)
jobs = [(f"act{i:02d}", 1, s) for i, s in enumerate(acts)] + \
       [(f"dec{i:03d}", 0, s) for i, s in enumerate(decs)]
random.shuffle(jobs)

done = set()
if os.path.exists("dude_scores.csv"):
    done = {r["id"] for r in csv.DictReader(open("dude_scores.csv"))}
newfile = not os.path.exists("dude_scores.csv")
fh = open("dude_scores.csv", "a", newline="")
wr = csv.writer(fh)
if newfile: wr.writerow(["id", "label", "smiles", "vina_score"])

for name, label, smi in jobs:
    if name in done: continue
    try:
        mol = Chem.AddHs(Chem.MolFromSmiles(smi))
        if AllChem.EmbedMolecule(mol, randomSeed=7) != 0: raise RuntimeError("embed")
        AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
        sdf = f"lig/{name}.sdf"
        w = Chem.SDWriter(sdf); w.write(mol); w.close()
        subprocess.run(["obabel", sdf, "-O", f"lig/{name}.pdbqt"], check=True, capture_output=True)
        out = f"out/{name}.pdbqt"
        subprocess.run(["vina", "--receptor", REC, "--ligand", f"lig/{name}.pdbqt",
            "--center_x", str(CENTER[0]), "--center_y", str(CENTER[1]), "--center_z", str(CENTER[2]),
            "--size_x", "22", "--size_y", "22", "--size_z", "22",
            "--exhaustiveness", "8", "--num_modes", "1", "--seed", "1", "--cpu", "2",
            "--out", out], check=True, capture_output=True, timeout=900)
        score = None
        for l in open(out):
            if l.startswith("REMARK VINA RESULT"):
                score = float(l.split()[3]); break
        wr.writerow([name, label, smi, score]); fh.flush()
        print(name, label, score, flush=True)
    except Exception as e:
        print("SKIPPED", name, label, type(e).__name__, flush=True)
