import os, sys, glob, subprocess, csv

CENTER = (22.01, 0.25, 52.79)   # erlotinib centre in 1M17
SIZE = 22
os.makedirs("results", exist_ok=True)

recs = sys.argv[1].split(",") if len(sys.argv) > 1 else \
    [os.path.basename(p)[:-6] for p in sorted(glob.glob("receptors/*.pdbqt"))]
ligs = sys.argv[2].split(",") if len(sys.argv) > 2 else \
    [os.path.basename(p)[:-6] for p in sorted(glob.glob("ligands/*.pdbqt"))]

def best_score(path):
    with open(path) as fh:
        for line in fh:
            if line.startswith("REMARK VINA RESULT"):
                return float(line.split()[3])

rows = []
for rec in recs:
    for lig in ligs:
        out = f"results/{rec}__{lig}.pdbqt"
        if not os.path.exists(out):
            subprocess.run(["vina", "--receptor", f"receptors/{rec}.pdbqt",
                "--ligand", f"ligands/{lig}.pdbqt",
                "--center_x", str(CENTER[0]), "--center_y", str(CENTER[1]),
                "--center_z", str(CENTER[2]),
                "--size_x", str(SIZE), "--size_y", str(SIZE), "--size_z", str(SIZE),
                "--exhaustiveness", "16", "--seed", "1", "--cpu", "2",
                "--out", out], check=True, capture_output=True)
        s = best_score(out)
        print(f"{rec:18s} {lig:12s} {s}", flush=True)
        rows.append((rec, lig, s))

with open("scores.csv", "a" if os.path.exists("scores.csv") else "w", newline="") as fh:
    wr = csv.writer(fh)
    if fh.tell() == 0:
        wr.writerow(["receptor", "ligand", "vina_score"])
    wr.writerows(rows)
