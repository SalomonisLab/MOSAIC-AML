#!/usr/bin/env python3
"""Calibration (reliability + Brier/ECE) and decision-curve analysis from the per-sample OOF scores.

Reads scratchpad/oof_scores_v2.json (multimodal arm: per-mutation score,y,donor). Runs locally.
  - Reliability: pooled (score,y) across mutations; raw vs nested-isotonic-calibrated; Brier + ECE.
  - Decision curve: per key mutation, net benefit vs threshold-probability, vs treat-all / treat-none.
Writes deliverables/calibration_reliability.{png,pdf} + decision_curves.{png,pdf} + calibration_dca.json.
"""
import os, json
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.isotonic import IsotonicRegression
from sklearn.model_selection import GroupKFold

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTD = os.path.join(ROOT, "deliverables"); os.makedirs(OUTD, exist_ok=True)
S = json.load(open(os.path.join(ROOT, "scratchpad", "oof_scores_v2.json")))["multimodal"]
BLUE, GREEN, RED, INK, MUTED, GRID, SURF = "#2a78d6", "#008300", "#c0392b", "#1a1a19", "#6b6b63", "#d9d9d4", "#ffffff"

# dedup aliases
def canon(m):
    ml=m.lower()
    return "inv16" if ("inv(16)" in ml or ml=="inv16") else ("kmt2a" if "kmt2a" in ml else m)
seen, MUTS = set(), []
for m in S:
    c = canon(m)
    if c not in seen: seen.add(c); MUTS.append(m)

# ---------- pooled reliability + Brier/ECE, raw vs nested-isotonic ----------
allp, ally, cal = [], [], []
CALM = {}                                                 # per-mutation nested-calibrated probability (for DCA)
for m in MUTS:
    p = np.array(S[m]["score"], float); y = np.array(S[m]["y"], int)
    g = np.array(S[m]["donor"])
    allp.append(p); ally.append(y)
    # nested isotonic: fit on other donor folds, apply to held fold (honest calibrated prob)
    cc = np.full(len(y), np.nan)
    ng = min(3, len(set(g)))
    if ng >= 2 and len(set(y)) == 2:
        for tri, tei in GroupKFold(ng).split(p, y, g):
            if len(set(y[tri])) < 2:
                cc[tei] = y[tri].mean(); continue
            iso = IsotonicRegression(out_of_bounds="clip").fit(p[tri], y[tri])
            cc[tei] = iso.predict(p[tei])
    else:
        cc = p.copy()
    CALM[m] = cc
    cal.append(cc)
allp = np.concatenate(allp); ally = np.concatenate(ally); cal = np.concatenate(cal)
ok = ~np.isnan(cal); allp, ally, cal = allp[ok], ally[ok], cal[ok]

def reliability(p, y, nb=10):
    edges = np.linspace(0, 1, nb + 1); xs, ys, ws = [], [], []
    for i in range(nb):
        m = (p >= edges[i]) & (p < edges[i + 1] if i < nb - 1 else p <= 1.0001)
        if m.sum() >= 5:
            xs.append(p[m].mean()); ys.append(y[m].mean()); ws.append(m.mean())
    return np.array(xs), np.array(ys), np.array(ws)

def brier(p, y): return float(np.mean((p - y) ** 2))
def ece(p, y, nb=10):
    edges = np.linspace(0, 1, nb + 1); e = 0.0
    for i in range(nb):
        m = (p >= edges[i]) & (p < edges[i + 1] if i < nb - 1 else p <= 1.0001)
        if m.sum(): e += m.mean() * abs(y[m].mean() - p[m].mean())
    return float(e)

rx, ry, rw = reliability(allp, ally); cx, cy, cw = reliability(cal, ally)
metrics = {"raw_brier": round(brier(allp, ally), 4), "cal_brier": round(brier(cal, ally), 4),
           "raw_ece": round(ece(allp, ally), 4), "cal_ece": round(ece(cal, ally), 4),
           "n": int(len(ally)), "prevalence": round(float(ally.mean()), 4)}

fig, ax = plt.subplots(1, 2, figsize=(12.5, 5.4), facecolor=SURF, gridspec_kw={"width_ratios": [1.35, 1]})
ax[0].plot([0, 1], [0, 1], color=MUTED, ls=(0, (5, 4)), lw=1.2, label="perfect")
ax[0].plot(rx, ry, "-o", color="#b9b9b0", lw=2, ms=6, label="raw fused score")
ax[0].plot(cx, cy, "-o", color=GREEN, lw=2.4, ms=6, label="isotonic-calibrated (nested)")
ax[0].set_xlabel("predicted probability", fontsize=10); ax[0].set_ylabel("observed frequency of mutation", fontsize=10)
ax[0].set_xlim(0, 1); ax[0].set_ylim(0, 1); ax[0].legend(frameon=False, fontsize=9, loc="upper left", labelcolor=INK)
ax[0].set_title("Reliability — pooled over %d mutations (n=%d)" % (len(MUTS), metrics["n"]), fontsize=12, color=INK, fontweight="bold")
for a in ax:
    a.set_facecolor(SURF); a.spines[["top","right"]].set_visible(False); a.spines[["left","bottom"]].set_color(GRID)
    a.tick_params(colors=MUTED, labelsize=9); a.grid(True, color=GRID, lw=.6); a.set_axisbelow(True)
xb = np.arange(2)
ax[1].bar(xb-0.19, [metrics["raw_brier"], metrics["raw_ece"]], 0.38, color="#b9b9b0", label="raw")
ax[1].bar(xb+0.19, [metrics["cal_brier"], metrics["cal_ece"]], 0.38, color=GREEN, label="calibrated")
for i,(r,c) in enumerate([(metrics["raw_brier"],metrics["cal_brier"]),(metrics["raw_ece"],metrics["cal_ece"])]):
    ax[1].text(i-0.19, r+0.003, "%.3f"%r, ha="center", va="bottom", fontsize=8.5, color=MUTED)
    ax[1].text(i+0.19, c+0.003, "%.3f"%c, ha="center", va="bottom", fontsize=8.5, color=GREEN, fontweight="bold")
ax[1].set_xticks(xb); ax[1].set_xticklabels(["Brier ↓","ECE ↓"], fontsize=10); ax[1].legend(frameon=False, fontsize=9, labelcolor=INK)
ax[1].set_title("Lower is better", fontsize=12, color=INK, fontweight="bold")
fig.suptitle("Prediction calibration — a stated probability matches observed frequency after calibration",
             fontsize=13.5, color=INK, fontweight="bold", y=0.99)
fig.subplots_adjust(left=0.07, right=0.98, top=0.86, bottom=0.11, wspace=0.2)
fig.savefig(os.path.join(OUTD, "calibration_reliability.pdf"), facecolor=SURF)
fig.savefig(os.path.join(OUTD, "calibration_reliability.png"), dpi=130, facecolor=SURF); plt.close(fig)

# ---------- decision-curve analysis, per key mutation ----------
KEY = [m for m in ["NPM1","inv16","TP53","TET2","FLT3-ITD","complex"] if m in {canon(x):x for x in MUTS} or m in S]
KEY = [next((x for x in MUTS if canon(x)==k or x==k), None) for k in ["NPM1","inv16","TP53","TET2","FLT3-ITD","complex"]]
KEY = [k for k in KEY if k]
pts = np.linspace(0.01, 0.6, 60)
dca = {}
fig, axes = plt.subplots(2, 3, figsize=(14, 7.6), facecolor=SURF); axes = axes.ravel()
for ax_i, m in enumerate(KEY):
    prob = CALM[m]; y = np.array(S[m]["y"], int); n = len(y); prev = y.mean()
    ok = ~np.isnan(prob); prob, yv = prob[ok], y[ok]; n = len(yv); prev = yv.mean()
    nb_model, nb_all = [], []
    for pt in pts:
        call = prob >= pt                                 # calibrated probability thresholded at pt (correct DCA)
        tp = np.sum(call & (yv == 1)); fp = np.sum(call & (yv == 0))
        nb_model.append(tp/n - (fp/n)*(pt/(1-pt)))
        nb_all.append(prev - (1-prev)*(pt/(1-pt)))
    ax = axes[ax_i]
    ax.plot(pts, nb_model, color=GREEN, lw=2.4, label="model")
    ax.plot(pts, nb_all, color=MUTED, lw=1.3, ls=(0,(4,3)), label="treat all")
    ax.axhline(0, color=RED, lw=1.1, ls=(0,(4,3)), label="treat none")
    ax.set_title("%s  (prev %.0f%%)" % (m, prev*100), fontsize=11, color=INK, fontweight="bold")
    ax.set_xlim(0, 0.6); ax.set_ylim(min(-0.02, min(nb_model)-0.01), max(0.02, prev*1.1))
    ax.set_facecolor(SURF); ax.spines[["top","right"]].set_visible(False); ax.spines[["left","bottom"]].set_color(GRID)
    ax.tick_params(colors=MUTED, labelsize=8.5); ax.grid(True, color=GRID, lw=.6); ax.set_axisbelow(True)
    if ax_i == 0: ax.legend(frameon=False, fontsize=8.5, loc="upper right", labelcolor=INK)
    if ax_i % 3 == 0: ax.set_ylabel("net benefit", fontsize=9.5)
    if ax_i >= 3: ax.set_xlabel("threshold probability", fontsize=9.5)
    dca[m] = {"pt": [round(x,3) for x in pts], "net_benefit": [round(float(x),4) for x in nb_model]}
fig.suptitle("Decision-curve analysis — where a call beats 'treat all' and 'treat none'", fontsize=14, color=INK, fontweight="bold", y=0.99)
fig.text(0.5, 0.01, "Net benefit vs threshold probability, per mutation. The model curve above both references = a positive call adds clinical value in that threshold range.",
         ha="center", fontsize=8.5, color=MUTED, style="italic")
fig.subplots_adjust(left=0.06, right=0.98, top=0.90, bottom=0.09, hspace=0.32, wspace=0.2)
fig.savefig(os.path.join(OUTD, "decision_curves.pdf"), facecolor=SURF)
fig.savefig(os.path.join(OUTD, "decision_curves.png"), dpi=130, facecolor=SURF); plt.close(fig)

json.dump({"calibration": metrics, "dca_mutations": list(dca)}, open(os.path.join(OUTD, "calibration_dca.json"), "w"), indent=1)
print("calibration:", metrics)
print("Brier %.3f -> %.3f  ECE %.3f -> %.3f (raw -> calibrated)" % (metrics["raw_brier"], metrics["cal_brier"], metrics["raw_ece"], metrics["cal_ece"]))
print("wrote calibration_reliability.* + decision_curves.* + calibration_dca.json")
print("CALIBRATION DCA OK")
