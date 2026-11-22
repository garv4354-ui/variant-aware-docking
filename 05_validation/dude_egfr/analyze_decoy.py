import pandas as pd, numpy as np, os
from sklearn.metrics import roc_auc_score, roc_curve
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

df = pd.read_csv("dude_scores.csv").dropna(subset=["vina_score"])
y = df["label"].values; s = -df["vina_score"].values     # higher = better binder
auc = roc_auc_score(y, s)
rng = np.random.default_rng(1); boots = []
for _ in range(2000):
    idx = rng.integers(0, len(y), len(y))
    if 0 < y[idx].sum() < len(idx):
        boots.append(roc_auc_score(y[idx], s[idx]))
lo, hi = np.percentile(boots, [2.5, 97.5])
top = max(1, round(0.1 * len(y)))
ef10 = y[np.argsort(-s)][:top].mean() / y.mean()
print(f"docked: {len(y)} (actives {int(y.sum())}, decoys {int((1-y).sum())})")
print(f"mean Vina score  actives {df[df.label==1].vina_score.mean():.2f}   decoys {df[df.label==0].vina_score.mean():.2f}")
print(f"ROC-AUC = {auc:.2f}  (95% bootstrap CI {lo:.2f} to {hi:.2f})")
print(f"Enrichment in top 10% = {ef10:.2f}")
fpr, tpr, _ = roc_curve(y, s)
os.makedirs("../../results/05_validation", exist_ok=True)
plt.figure(figsize=(5, 5)); plt.plot(fpr, tpr, label=f"AUC = {auc:.2f}")
plt.plot([0, 1], [0, 1], "k--", label="random"); plt.xlabel("False positive rate")
plt.ylabel("True positive rate"); plt.title("DUD-E EGFR subset, Vina on EGFR WT (1M17)")
plt.legend(); plt.tight_layout(); plt.savefig("../../results/05_validation/roc_dude.png", dpi=200)
