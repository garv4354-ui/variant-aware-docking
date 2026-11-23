import os, csv, subprocess, statistics as st
from pdbfixer import PDBFixer
from openmm.app import PDBFile

CENTER = (22.01, 0.25, 52.79)
SIZE = 22
LIGS = ["erlotinib", "osimertinib", "aspirin"]
SEEDS = [1, 2, 3, 4, 5]
os.makedirs("results_af", exist_ok=True)

# 1. Prepare the AlphaFold receptor the same way as the experimental ones
fixer = PDBFixer(filename="../03_structures/AF_aligned_trimmed.pdb")
fixer.removeHeterogens(keepWater=False)
fixer.findMissingResidues(); fixer.missingResidues = {}
fixer.findNonstandardResidues(); fixer.replaceNonstandardResidues()
fixer.findMissingAtoms(); fixer.addMissingAtoms()
fixer.addMissingHydrogens(7.4)
with open("../03_structures/AF_prepared.pdb", "w") as fh:
    PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)
subprocess.run(["obabel", "-ipdb", "../03_structures/AF_prepared.pdb",
                "-xr", "-O", "receptors/EGFR_AF.pdbqt"],
               check=True, capture_output=True)
print("receptor ready: EGFR_AF", flush=True)

def best(path):
    with open(path) as fh:
        for line in fh:
            if line.startswith("REMARK VINA RESULT"):
                return float(line.split()[3])

# 2. Dock with the same settings as the seed test
rows, af_means = [], {}
for lig in LIGS:
    scores = []
    for seed in SEEDS:
        out = f"results_af/EGFR_AF__{lig}__s{seed}.pdbqt"
        if not os.path.exists(out):
            subprocess.run(["vina", "--receptor", "receptors/EGFR_AF.pdbqt",
                            "--ligand", f"ligands/{lig}.pdbqt",
                            "--center_x", str(CENTER[0]), "--center_y", str(CENTER[1]),
                            "--center_z", str(CENTER[2]),
                            "--size_x", str(SIZE), "--size_y", str(SIZE), "--size_z", str(SIZE),
                            "--exhaustiveness", "16", "--seed", str(seed), "--cpu", "2",
                            "--out", out], check=True, capture_output=True)
        scores.append(best(out))
    m, sd = st.mean(scores), st.stdev(scores)
    af_means[lig] = m
    rows.append(["EGFR_AF", lig, len(scores), round(m, 3), round(sd, 3), min(scores), max(scores)])
    print(f"EGFR_AF     {lig:12s} n={len(scores)} mean={m:.2f} sd={sd:.2f} "
          f"range={min(scores):.2f}..{max(scores):.2f}", flush=True)

# 3. Compare with experimental wild-type means from the seed test
wt = {}
with open("seed_scores.csv") as fh:
    for r in csv.DictReader(fh):
        if r["receptor"] == "EGFR_WT":
            wt[r["ligand"]] = float(r["mean"])
print("\nAlphaFold minus experimental WT (mean score):")
for lig in LIGS:
    print(f"  {lig:12s} {af_means[lig] - wt[lig]:+.2f} kcal/mol")

with open("af_scores.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["receptor", "ligand", "n", "mean", "sd", "min", "max"])
    w.writerows(rows)
print("DONE")
