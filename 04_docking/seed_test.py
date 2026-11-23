import os, csv, subprocess, statistics as st

CENTER = (22.01, 0.25, 52.79)
SIZE = 22
RECS = ["EGFR_WT", "EGFR_T790M"]
LIGS = ["erlotinib", "osimertinib", "aspirin"]
SEEDS = [2, 3, 4, 5]
os.makedirs("results_seed", exist_ok=True)

def best(path):
    with open(path) as fh:
        for line in fh:
            if line.startswith("REMARK VINA RESULT"):
                return float(line.split()[3])

means = {}
rows = []
for rec in RECS:
    for lig in LIGS:
        scores = []
        first = f"results/{rec}__{lig}.pdbqt"
        if os.path.exists(first):
            scores.append(best(first))
        for seed in SEEDS:
            out = f"results_seed/{rec}__{lig}__s{seed}.pdbqt"
            if not os.path.exists(out):
                subprocess.run(["vina", "--receptor", f"receptors/{rec}.pdbqt",
                    "--ligand", f"ligands/{lig}.pdbqt",
                    "--center_x", str(CENTER[0]), "--center_y", str(CENTER[1]),
                    "--center_z", str(CENTER[2]),
                    "--size_x", str(SIZE), "--size_y", str(SIZE), "--size_z", str(SIZE),
                    "--exhaustiveness", "16", "--seed", str(seed), "--cpu", "2",
                    "--out", out], check=True, capture_output=True)
            scores.append(best(out))
        m, sd = st.mean(scores), st.stdev(scores)
        means[(rec, lig)] = m
        rows.append([rec, lig, len(scores), round(m, 3), round(sd, 3),
                     min(scores), max(scores)])
        print(f"{rec:11s} {lig:12s} n={len(scores)} mean={m:.2f} sd={sd:.2f} "
              f"range={min(scores):.2f}..{max(scores):.2f}", flush=True)

print("\nMutant minus wild-type (mean score):")
for lig in LIGS:
    d = means[("EGFR_T790M", lig)] - means[("EGFR_WT", lig)]
    print(f"  {lig:12s} {d:+.2f} kcal/mol")

with open("seed_scores.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["receptor", "ligand", "n", "mean", "sd", "min", "max"])
    w.writerows(rows)
print("DONE")
