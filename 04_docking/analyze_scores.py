import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

df = pd.read_csv("scores.csv").drop_duplicates(["receptor", "ligand"], keep="last")
df["receptor"] = df["receptor"].str.replace("EGFR_", "", regex=False)
order = ["WT", "L858R", "T790M", "L858R_T790M", "C797S", "G719S"]
egfr = ["erlotinib", "gefitinib", "afatinib", "dacomitinib", "osimertinib", "lapatinib", "brigatinib"]
ctrl = ["aspirin", "caffeine", "ibuprofen"]
tab = df.pivot(index="ligand", columns="receptor", values="vina_score")[order].loc[egfr + ctrl]
delta = tab.drop(columns="WT").sub(tab["WT"], axis=0)   # positive = weaker in mutant

tab.round(2).to_csv("../results/04_docking/score_table.csv")
delta.round(2).to_csv("../results/04_docking/delta_vs_WT.csv")
print("Vina scores (kcal/mol):\n", tab.round(2), "\n")
print("Mutant minus WT (positive = weaker in mutant):\n", delta.round(2))

fig, ax = plt.subplots(figsize=(8, 6))
im = ax.imshow(delta.values, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
ax.set_xticks(range(delta.shape[1])); ax.set_xticklabels(delta.columns, rotation=30, ha="right")
ax.set_yticks(range(delta.shape[0])); ax.set_yticklabels(delta.index)
for i in range(delta.shape[0]):
    for j in range(delta.shape[1]):
        ax.text(j, i, f"{delta.values[i, j]:+.2f}", ha="center", va="center", fontsize=8)
ax.axhline(len(egfr) - 0.5, color="black", lw=1.5)
ax.set_title("Vina score change vs wild-type (kcal/mol)\nred = weaker binding, blue = stronger")
fig.colorbar(im, ax=ax); fig.tight_layout()
fig.savefig("../results/04_docking/heatmap_delta.png", dpi=200)
print("saved heatmap")
