#!/usr/bin/env python3
"""Extract the full per-(mutation x modality) CV-OOF AUROC breakdown from the DEPLOYED predictor
(mutation_predictor.pkl) — the optimized sc system (StandardScaler->diff_select(500)->LinearSVC(C=0.02)
->percentile, donor-grouped CV-OOF, ridge-NNLS fusion, F1-max thresholds). No retrain."""
import os, sys, json
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
ROOT = r"C:\Users\krog5w\.gemini\antigravity\scratch\AML-multimodal"
sys.path.insert(0, os.path.join(ROOT, "pipeline"))
from amlmm.predictor import MutationPredictor

P = MutationPredictor.load(os.path.join(ROOT, "pipeline", "mutation_predictor.pkl"))
card = json.load(open(os.path.join(ROOT, "pipeline", "model_card.json")))["mutations"]
MODS = P.modalities
print("deployed: %d mutations x %d modalities: %s" % (len(P.mutations), len(MODS), ", ".join(MODS)))

rows = []
for m in P.mutations:
    r = {"mutation": m, "n_pos": card.get(m, {}).get("n_pos"),
         "fused_train_auc": P.train_auc.get(m), "heldout_auc": P.heldout_auc.get(m),
         "threshold": P.thresholds.get(m)}
    percv = {}
    for mod in MODS:
        mm = P.models.get((m, mod))
        r[mod] = mm.oof_auc if (mm is not None and mm.oof_auc is not None) else np.nan
        if mm is not None and mm.oof_auc is not None:
            percv[mod] = mm.oof_auc
    r["best_modality"] = max(percv, key=percv.get) if percv else None
    r["best_single_auc"] = max(percv.values()) if percv else np.nan
    w = P.weights.get(m, {})
    r["fusion_weights"] = {k: v for k, v in sorted(w.items(), key=lambda x: -x[1]) if v > 0}
    rows.append(r)
df = pd.DataFrame(rows).sort_values("fused_train_auc", ascending=False).reset_index(drop=True)
df.to_csv(os.path.join(ROOT, "scratchpad", "deployed_modality_breakdown.csv"), index=False)

# ---- which modality carries the system ----
print("\n== per-modality: mean CV-OOF AUROC over mutations | #times best | #times weighted>0 ==")
for mod in MODS:
    col = df[mod].dropna()
    nbest = int((df["best_modality"] == mod).sum())
    nwt = int(sum(1 for w in df["fusion_weights"] if mod in w))
    print("  %-11s mean=%.3f (n=%2d)  best=%2d  weighted=%2d" % (mod, col.mean(), len(col), nbest, nwt))

# ---- fusion vs best single: does the NNLS combiner beat the best modality? ----
d2 = df.dropna(subset=["fused_train_auc", "best_single_auc"])
gain = (d2["fused_train_auc"] - d2["best_single_auc"])
print("\n== fusion vs best-single-modality (train CV-OOF) ==")
print("  mean fused=%.3f  mean best-single=%.3f  mean gain=%+.3f  fusion>=best in %d/%d muts"
      % (d2["fused_train_auc"].mean(), d2["best_single_auc"].mean(), gain.mean(), int((gain >= -1e-9).sum()), len(d2)))

# ---- FLT3-ITD (the one variant-level split feasible in sc) ----
for key in ["FLT3-ITD", "FLT3_ITD", "FLT3"]:
    if key in df["mutation"].values:
        rr = df[df["mutation"] == key].iloc[0]
        print("\n== %s per-modality CV-OOF AUROC ==" % key)
        print("   " + "  ".join("%s=%.2f" % (mod, rr[mod]) for mod in MODS if rr[mod] == rr[mod]))
        print("   fused=%.3f heldout=%s best=%s" % (rr["fused_train_auc"], rr["heldout_auc"], rr["best_modality"]))
        break

# =========================== FIGURE ===========================
plt.rcParams.update({"font.size": 8.5, "figure.dpi": 130})
fig = plt.figure(figsize=(15, 9.5))
gs = fig.add_gridspec(1, 2, width_ratios=[1.5, 1], wspace=0.24)
# Panel A: heatmap mutation x modality
axA = fig.add_subplot(gs[0, 0])
Mx = df[MODS].values.astype(float)
im = axA.imshow(Mx, aspect="auto", cmap="RdYlGn", vmin=0.5, vmax=1.0)
axA.set_xticks(range(len(MODS))); axA.set_xticklabels(MODS, rotation=35, ha="right")
axA.set_yticks(range(len(df))); axA.set_yticklabels(["%s (n%s)" % (m, n) for m, n in zip(df["mutation"], df["n_pos"])], fontsize=7)
for i in range(len(df)):
    for j, mod in enumerate(MODS):
        v = Mx[i, j]
        if v == v:
            axA.text(j, i, "%.2f" % v, ha="center", va="center", fontsize=6.2,
                     color="black" if 0.6 < v < 0.9 else "white")
        if df.iloc[i]["best_modality"] == mod:                 # box the winning modality
            axA.add_patch(Rectangle((j - .5, i - .5), 1, 1, fill=False, edgecolor="#111", lw=1.8))
    # star modalities that get fusion weight
    for mod in df.iloc[i]["fusion_weights"]:
        j = MODS.index(mod); axA.text(j, i + 0.32, "*", ha="center", va="center", fontsize=8, color="#0033aa")
axA.set_title("A. Per-modality CV-OOF AUROC (deployed sc system)  |  □ = best modality,  * = used in NNLS fusion",
              fontsize=9.5, loc="left")
cb = fig.colorbar(im, ax=axA, fraction=0.025, pad=0.01); cb.set_label("CV-OOF AUROC")
# Panel B: fused vs best-single vs heldout
axB = fig.add_subplot(gs[0, 1])
y = np.arange(len(df))
axB.barh(y + 0.22, df["best_single_auc"], height=0.28, color="#95a5a6", label="best single modality")
axB.barh(y - 0.06, df["fused_train_auc"], height=0.28, color="#2980b9", label="fused (NNLS, train CV-OOF)")
axB.scatter(df["heldout_auc"], y - 0.06, s=16, color="#c0392b", zorder=5, label="sealed held-out")
axB.set_yticks(y); axB.set_yticklabels(df["mutation"], fontsize=7); axB.invert_yaxis()
axB.set_xlim(0.5, 1.03); axB.axvline(0.5, color="#ccc", ls=":")
axB.set_xlabel("AUROC"); axB.legend(fontsize=7, loc="lower right")
axB.set_title("B. Fusion vs best single modality\n+ sealed held-out (●)", fontsize=9.5, loc="left")
fig.suptitle("Deployed multimodal sc mutation predictor — per-modality breakdown  |  26 mutations x 8 sc-derived modalities  |  all optimizations (diff-select 500, LinearSVC C=0.02, donor-grouped CV-OOF, ridge-NNLS fusion, F1-max thresholds)",
             fontsize=10.5, y=0.995)
out1 = os.path.join(ROOT, "scratchpad", "deployed_modality_breakdown.png")
fig.savefig(out1, bbox_inches="tight")
print("\nwrote", out1); print("wrote deployed_modality_breakdown.csv")
