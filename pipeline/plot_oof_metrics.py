#!/usr/bin/env python3
"""PDF deliverables: multimodal vs bulk RNA-seq alone (no imputation), sensitivity / specificity / F1.

Reads scratchpad/oof_metrics.json (from eval_oof_metrics.py) and writes, into deliverables/:
  1. ..._overview.pdf      all mutations together — mean bar per arm with ONE DOT PER MUTATION
                           overlaid and paired lines connecting each mutation's two arms.
  2. ..._per_mutation.pdf  one panel per mutation, grouped bars (sens / spec / F1) x 2 arms.
  3. ..._metrics.tsv       the underlying numbers.

Both arms are donor-grouped CV-OOF on the SAME samples, labels and folds — differing only in which
modality blocks the model was given, so the contrast is paired and fair.

  /usr/local/anaconda3-2020/bin/python plot_oof_metrics.py
"""
import os, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SRC = os.path.join(ROOT, "scratchpad", "oof_metrics.json")
OUTD = os.path.join(ROOT, "deliverables")
os.makedirs(OUTD, exist_ok=True)

# validated categorical slots 1 & 2 (light surface) — see dataviz validator: CVD dE 26.5, normal 29.0
C_BULK, C_MULTI = "#2a78d6", "#008300"
SURFACE, INK, MUTED, GRID = "#ffffff", "#1a1a19", "#6b6b63", "#d9d9d4"
ARM_LABEL = {"bulkrna": "Bulk RNA-seq alone (no imputation)",
             "multimodal": "MOSAIC-AML multimodal (8 modalities)"}
METRICS = [("sensitivity", "Sensitivity"), ("specificity", "Specificity"), ("f1", "F1")]

d = json.load(open(SRC))
A = d["arms"]; S = d["paired_summary"]
shared = S["shared"]
gen = d.get("generated", "")


def _sig(m):
    """metric fingerprint across both arms — alias rows (inv16/inv(16)_CBFB-MYH11, kmt2a/
    KMT2A-rearrangement) are the SAME event scored twice and would double-count in the means."""
    return tuple(A[arm]["mutations"][m].get(k) for arm in ("bulkrna", "multimodal")
                 for k in ("sensitivity", "specificity", "f1", "auroc", "n_pos"))


_seen, _drop = {}, []
for m in sorted(shared):
    s = _sig(m)
    if s in _seen:                       # identical in every metric -> alias; keep the descriptive name
        keep, dup = (m, _seen[s]) if len(m) > len(_seen[s]) else (_seen[s], m)
        _seen[s] = keep; _drop.append(dup)
    else:
        _seen[s] = m
shared = [m for m in shared if m not in set(_drop)]
if _drop:
    print("collapsed %d alias row(s): %s" % (len(_drop), ", ".join(sorted(_drop))))

# recompute the paired summary on the DEduplicated set (the JSON's was over all rows incl. aliases)
S = {}
for key, _lab in [("sensitivity", ""), ("specificity", ""), ("f1", ""), ("auroc", "")]:
    b = [A["bulkrna"]["mutations"][m][key] for m in shared]
    a = [A["multimodal"]["mutations"][m][key] for m in shared]
    b = [x for x in b if x is not None]; a = [x for x in a if x is not None]
    S[key] = {"bulkrna_mean": float(np.mean(b)), "multimodal_mean": float(np.mean(a)),
              "delta": float(np.mean(a) - np.mean(b))}
    try:
        from scipy.stats import wilcoxon
        S[key]["wilcoxon_p"] = float(wilcoxon(
            [A["multimodal"]["mutations"][m][key] for m in shared],
            [A["bulkrna"]["mutations"][m][key] for m in shared])[1])
    except Exception:
        pass
CAVEAT = ("Donor-grouped 3-fold cross-validated out-of-fold predictions; both arms use identical samples, "
          "labels and folds.\nOperating point = F1-max threshold on the same OOF vector, so absolute values "
          "are optimistic — the CONTRAST between arms is the result.")


def style(ax):
    ax.set_facecolor(SURFACE)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color(GRID); ax.spines["bottom"].set_color(GRID)
    ax.tick_params(colors=MUTED, labelsize=9, length=3)
    ax.yaxis.grid(True, color=GRID, lw=0.7, alpha=0.8); ax.set_axisbelow(True)
    ax.xaxis.grid(False)


def val(arm, mut, key):
    v = A[arm]["mutations"][mut].get(key)
    return np.nan if v is None else float(v)


# ------------------------------------------------------------------ 1. overview
fig, axes = plt.subplots(1, 3, figsize=(13.5, 5.6), facecolor=SURFACE)
rng = np.random.default_rng(7)                                  # deterministic jitter
for ax, (key, title) in zip(axes, METRICS):
    b = np.array([val("bulkrna", m, key) for m in shared])
    a = np.array([val("multimodal", m, key) for m in shared])
    mb, ma = np.nanmean(b), np.nanmean(a)
    # mean bars — 2px surface gap via white edge
    ax.bar([0], [mb], width=0.52, color=C_BULK, alpha=0.20, edgecolor=SURFACE, lw=2, zorder=1)
    ax.bar([1], [ma], width=0.52, color=C_MULTI, alpha=0.20, edgecolor=SURFACE, lw=2, zorder=1)
    ax.plot([-0.26, 0.26], [mb, mb], color=C_BULK, lw=2.5, zorder=4, solid_capstyle="round")
    ax.plot([0.74, 1.26], [ma, ma], color=C_MULTI, lw=2.5, zorder=4, solid_capstyle="round")
    # paired dots — one per mutation, connected so the per-mutation direction is visible
    jb = rng.uniform(-0.10, 0.10, len(shared)); ja = rng.uniform(-0.10, 0.10, len(shared))
    for i in range(len(shared)):
        if np.isnan(b[i]) or np.isnan(a[i]):
            continue
        ax.plot([0 + jb[i], 1 + ja[i]], [b[i], a[i]], color=MUTED, lw=0.6, alpha=0.30, zorder=2)
    ax.scatter(0 + jb, b, s=46, color=C_BULK, edgecolor=SURFACE, lw=1.2, zorder=3, alpha=0.95)
    ax.scatter(1 + ja, a, s=46, color=C_MULTI, edgecolor=SURFACE, lw=1.2, zorder=3, alpha=0.95)
    # selective direct labels: the two means only (never a number on every dot)
    ax.text(0, -0.085, "%.3f" % mb, ha="center", va="top", fontsize=11, fontweight="bold", color=C_BULK)
    ax.text(1, -0.085, "%.3f" % ma, ha="center", va="top", fontsize=11, fontweight="bold", color=C_MULTI)
    p = S.get(key, {}).get("wilcoxon_p")
    sub = "Δ %+.3f" % (ma - mb) + ("   Wilcoxon p = %.2g" % p if p is not None else "")
    ax.set_title(title, fontsize=13, color=INK, fontweight="bold", pad=14)
    ax.text(0.5, 1.015, sub, transform=ax.transAxes, ha="center", va="bottom",
            fontsize=9.5, color=MUTED)
    ax.set_xticks([0, 1]); ax.set_xticklabels([])
    ax.set_xlim(-0.55, 1.55); ax.set_ylim(0, 1.04)
    style(ax)
axes[0].set_ylabel("score  (dot = one mutation)", fontsize=10, color=INK)

handles = [Patch(facecolor=C_BULK, edgecolor="none", label=ARM_LABEL["bulkrna"]),
           Patch(facecolor=C_MULTI, edgecolor="none", label=ARM_LABEL["multimodal"])]
# legend lives at the TOP so it can never collide with the mean labels under each bar
fig.legend(handles=handles, loc="upper center", ncol=2, frameon=False, fontsize=10.5,
           bbox_to_anchor=(0.5, 0.895), labelcolor=INK)
fig.suptitle("Multimodal vs bulk RNA-seq alone — %d driver mutations" % len(shared),
             fontsize=15.5, color=INK, fontweight="bold", y=0.982)
fig.text(0.5, 0.933, "paired donor-grouped CV-OOF on the single-cell AML atlas · identical samples & folds",
         ha="center", fontsize=10, color=MUTED)
fig.text(0.5, 0.022, CAVEAT, ha="center", fontsize=7.6, color=MUTED, style="italic")
fig.subplots_adjust(left=0.07, right=0.98, top=0.78, bottom=0.17, wspace=0.18)
p1 = os.path.join(OUTD, "MOSAIC-AML_multimodal_vs_bulkRNA_overview.pdf")
fig.savefig(p1, format="pdf", facecolor=SURFACE)
fig.savefig(p1.replace(".pdf", ".png"), format="png", dpi=130, facecolor=SURFACE)   # for visual QA
plt.close(fig)
print("wrote", p1)

# ------------------------------------------------------------------ 2. per mutation
order = sorted(shared, key=lambda m: -(val("multimodal", m, "f1") - val("bulkrna", m, "f1")))
ncol = 4
nrow = int(np.ceil(len(order) / ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(15, 2.65 * nrow), facecolor=SURFACE)
axes = np.array(axes).reshape(-1)
xs = np.arange(3); w = 0.36
for i, mut in enumerate(order):
    ax = axes[i]
    b = [val("bulkrna", mut, k) for k, _ in METRICS]
    a = [val("multimodal", mut, k) for k, _ in METRICS]
    ax.bar(xs - w / 2, b, w, color=C_BULK, edgecolor=SURFACE, lw=2, zorder=2)
    ax.bar(xs + w / 2, a, w, color=C_MULTI, edgecolor=SURFACE, lw=2, zorder=2)
    npos = A["multimodal"]["mutations"][mut].get("n_pos")
    ax.set_title("%s   (n+ = %s)" % (mut, npos), fontsize=10.5, color=INK, fontweight="bold", pad=6)
    ax.set_xticks(xs); ax.set_xticklabels([t for _, t in METRICS], fontsize=8.5)
    ax.set_ylim(0, 1.04); ax.set_yticks([0, 0.5, 1.0])
    style(ax)
    if i % ncol:
        ax.set_yticklabels([])
for j in range(len(order), len(axes)):
    axes[j].axis("off")
fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, fontsize=11,
           bbox_to_anchor=(0.5, 0.038), labelcolor=INK)
fig.suptitle("Per-mutation: sensitivity · specificity · F1", fontsize=15.5, color=INK,
             fontweight="bold", y=0.995)
fig.text(0.5, 0.977, "sorted by F1 gain (multimodal − bulk RNA-seq alone)",
         ha="center", fontsize=9.5, color=MUTED)
fig.text(0.5, 0.012, CAVEAT, ha="center", fontsize=8, color=MUTED, style="italic")   # wrapped, at foot
fig.subplots_adjust(left=0.055, right=0.985, top=0.952, bottom=0.075, hspace=0.62, wspace=0.12)
p2 = os.path.join(OUTD, "MOSAIC-AML_multimodal_vs_bulkRNA_per_mutation.pdf")
fig.savefig(p2, format="pdf", facecolor=SURFACE)
fig.savefig(p2.replace(".pdf", ".png"), format="png", dpi=110, facecolor=SURFACE)   # for visual QA
plt.close(fig)
print("wrote", p2)

# ------------------------------------------------------------------ 3. the numbers
p3 = os.path.join(OUTD, "MOSAIC-AML_multimodal_vs_bulkRNA_metrics.tsv")
with open(p3, "w") as fh:
    fh.write("mutation\tn\tn_pos\tprevalence\t" + "\t".join(
        "%s_%s" % (arm, k) for arm in ("bulkrna", "multimodal")
        for k in ("sensitivity", "specificity", "f1", "auroc", "auprc", "threshold")) + "\n")
    for mut in order:
        mm = A["multimodal"]["mutations"][mut]
        cells = [mut, str(mm["n"]), str(mm["n_pos"]), str(mm["prevalence"])]
        for arm in ("bulkrna", "multimodal"):
            r = A[arm]["mutations"][mut]
            cells += ["" if r.get(k) is None else str(r.get(k))
                      for k in ("sensitivity", "specificity", "f1", "auroc", "auprc", "threshold")]
        fh.write("\t".join(cells) + "\n")
    fh.write("\nMEAN over %d shared mutations\n" % len(shared))
    for k, lab in METRICS + [("auroc", "AUROC")]:
        s = S.get(k, {})
        fh.write("%s\tbulk %.4f\tmultimodal %.4f\tdelta %+.4f%s\n" %
                 (lab, s.get("bulkrna_mean", float("nan")), s.get("multimodal_mean", float("nan")),
                  s.get("delta", float("nan")),
                  ("\twilcoxon_p %.3g" % s["wilcoxon_p"]) if s.get("wilcoxon_p") is not None else ""))
print("wrote", p3)
print("PLOT OOF METRICS OK")
