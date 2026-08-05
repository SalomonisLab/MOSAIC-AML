#!/usr/bin/env python3
"""Per-mutation performance scatter panels — one dot per mutation, panels by (method x dataset).

Two companion views:
  sens_spec_panels.*        sensitivity (X) vs specificity (Y)   [as requested]
  precision_recall_panels.* recall (X) vs precision (Y)          [journal-figure style]

Specificity saturates near 1 on rare mutations (every model calls most samples absent correctly), so
the sens/spec panels are y-zoomed to the occupied band; precision-recall is the more discriminating
view for rare classes and is provided alongside.

  /usr/local/anaconda3-2020/bin/python build_sens_spec_panels.py    (runs locally too)
"""
import os, json, glob
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTD = os.path.join(ROOT, "deliverables"); os.makedirs(OUTD, exist_ok=True)
BLUE, GREEN, INK, MUTED, GRID, SURF, STAR = "#2a78d6", "#008300", "#1a1a19", "#6b6b63", "#d9d9d4", "#ffffff", "#b06a00"

# ---------- gather (method, dataset) -> {mut: dict(sens, spec, prec, n_pos)} ----------
panels = []
mm = json.load(open(os.path.join(ROOT, "scratchpad", "oof_metrics_v3_nyu2.json")))["arms"]["multimodal"]["mutations"]
seen, d = set(), {}
for m, r in mm.items():                                   # dedup alias rows (inv16 / inv(16)_CBFB-MYH11 etc.)
    key = (round(r["auroc"], 4), r["n_pos"])
    if key in seen: continue
    seen.add(key)
    tp, fp, fn = r["tp"], r["fp"], r["fn"]
    d[m] = dict(sens=r["nested_sensitivity"], spec=r["nested_specificity"],
                prec=(tp / (tp + fp) if tp + fp else 0.0), n=r["n_pos"])
panels.append(("MOSAIC-AML multimodal", "all scRNA (CV-OOF)", GREEN, d))

conf = {}                                                  # multimodal on the sealed held-out
for f in sorted(glob.glob(os.path.join(ROOT, "runs", "predict_*", "patient_report.json"))):
    rep = json.load(open(f))
    for p in (rep.get("mutation_predictions") or []):
        tl = p.get("true_label")
        if tl is None: continue
        c = conf.setdefault(p["mutation"], [0, 0, 0, 0])
        pos, cp = (tl == "present"), (p.get("call") == "present")
        c[0 if (pos and cp) else 2 if pos else 1 if cp else 3] += 1
d = {}
for m, (tp, fp, fn, tn) in conf.items():
    if tp + fn >= 1 and tn + fp >= 1:
        d[m] = dict(sens=tp / (tp + fn), spec=tn / (tn + fp),
                    prec=(tp / (tp + fp) if tp + fp else 0.0), n=tp + fn)
panels.append(("MOSAIC-AML multimodal", "held-out scRNA (n=29)", GREEN, d))

bm = json.load(open(os.path.join(OUTD, "bulk_matrix.json")))["cohorts"]
for coh, lab in [("BeatAML_CV", "BeatAML (CV, bulk)"), ("Leucegene", "Leucegene (external bulk)"),
                 ("all_scRNA", "all scRNA (cross-platform)")]:
    pc = bm[coh]["per_category"]
    d = {c: dict(sens=v["sensitivity"], spec=v["specificity"], prec=v.get("precision") or 0.0, n=v["n_pos"])
         for c, v in pc.items() if v.get("sensitivity") is not None and v.get("specificity") is not None}
    panels.append(("Bulk RNA-only model", lab, BLUE, d))


def draw(mode):
    """mode='spec' -> sensitivity vs specificity ; mode='prec' -> recall vs precision"""
    ykey, ylab = ("spec", "specificity") if mode == "spec" else ("prec", "precision")
    ncol = 3; nrow = int(np.ceil(len(panels) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(14.5, 4.8 * nrow), facecolor=SURF)
    axes = np.array(axes).reshape(-1)
    # shared y-range: zoom to the occupied band so differences are visible
    allv = [v[ykey] for _, _, _, dd in panels for v in dd.values()]
    ylo = max(0.0, min(allv) - 0.04) if mode == "spec" else 0.0
    for ax, (meth, ds, col, dd) in zip(axes, panels):
        xs = [v["sens"] for v in dd.values()]; ys = [v[ykey] for v in dd.values()]
        sz = [28 + 3.2 * v["n"] for v in dd.values()]
        ax.scatter(xs, ys, s=sz, color=col, edgecolor=SURF, lw=1.1, alpha=0.82, zorder=3)
        mx, my = float(np.mean(xs)), float(np.mean(ys))
        ax.axvline(mx, color=col, lw=1, ls=(0, (4, 3)), alpha=.5); ax.axhline(my, color=col, lw=1, ls=(0, (4, 3)), alpha=.5)
        if 1.0 <= 1.02 and ylo < 1.0:
            ax.plot([1], [1], marker="*", ms=13, color=STAR, zorder=4, clip_on=False)
        # label only the 3 worst + 2 best by distance to (1,1); skip the crowded top-right cluster
        order = sorted(dd, key=lambda m: ((dd[m]["sens"] - 1) ** 2 + (dd[m][ykey] - 1) ** 2))
        for m, off in [(order[0], (5, 5))] + [(m, (5, -9)) for m in order[-3:]]:
            ax.annotate(str(m)[:13], (dd[m]["sens"], dd[m][ykey]), fontsize=7.2, color=INK,
                        xytext=off, textcoords="offset points", zorder=5)
        ax.set_xlim(-0.03, 1.06); ax.set_ylim(ylo, 1.0 + (1.0 - ylo) * 0.06)
        ax.set_title("%s\n%s  ·  %d mutations" % (meth, ds, len(dd)), fontsize=10.5, color=INK, fontweight="bold", pad=8)
        ax.set_xlabel("sensitivity (recall)", fontsize=9.5); ax.set_ylabel(ylab, fontsize=9.5)
        ax.set_facecolor(SURF); ax.spines[["top", "right"]].set_visible(False)
        ax.spines[["left", "bottom"]].set_color(GRID); ax.tick_params(colors=MUTED, labelsize=8.5)
        ax.grid(True, color=GRID, lw=.6); ax.set_axisbelow(True)
        ax.text(0.02, ylo + (1 - ylo) * 0.06, "mean  sens %.2f\n         %s %.2f" % (mx, ylab[:4], my),
                fontsize=8.2, color=col, fontweight="bold", va="bottom")
    for j in range(len(panels), len(axes)):
        axes[j].axis("off")
    if len(panels) < len(axes):
        a = axes[len(panels)]; a.axis("off")
        for i, n in enumerate([5, 20, 50]):
            a.scatter([0.16], [0.68 - i * 0.13], s=28 + 3.2 * n, color=MUTED, alpha=.45, edgecolor=SURF)
            a.text(0.29, 0.68 - i * 0.13, "n⁺ = %d positives" % n, fontsize=9, va="center", color=INK)
        a.text(0.04, 0.88, "Each dot = one mutation", fontsize=11, fontweight="bold", color=INK)
        a.text(0.04, 0.26, ("Top-right = ideal.  Dot size = positives\navailable.  Dashed = panel means.\n"
                            + ("y-axis zoomed: specificity saturates\nnear 1 because drivers are rare."
                               if mode == "spec" else
                               "Precision is the discriminating axis\nfor rare mutations (few positives,\nmany chances to false-positive).")),
               fontsize=8.6, color=MUTED, va="top")
        a.set_xlim(0, 1); a.set_ylim(0, 1)
    ttl = ("Per-mutation sensitivity vs specificity, by model and dataset" if mode == "spec"
           else "Per-mutation precision vs recall, by model and dataset")
    fig.suptitle(ttl, fontsize=15, color=INK, fontweight="bold", y=0.985)
    fig.text(0.5, 0.008, "Multimodal at its nested-CV (honest) threshold; bulk model at its deployed F1-max threshold. "
             + ("Specificity is uniformly high — sensitivity separates the models."
                if mode == "spec" else
                "Precision falls with prevalence, so rare drivers sit low even when their ranking (AUROC) is strong."),
             ha="center", fontsize=8, color=MUTED, style="italic")
    fig.subplots_adjust(left=0.055, right=0.985, top=0.90, bottom=0.075, hspace=0.38, wspace=0.22)
    name = "sens_spec_panels" if mode == "spec" else "precision_recall_panels"
    for ext, dpi in [("pdf", None), ("png", 130)]:
        fig.savefig(os.path.join(OUTD, "%s.%s" % (name, ext)), **({"dpi": dpi} if dpi else {}), facecolor=SURF)
    plt.close(fig)
    print("wrote deliverables/%s.pdf/.png" % name)


draw("spec"); draw("prec")
print("panels:", [(m, ds, len(x)) for m, ds, _, x in panels])
