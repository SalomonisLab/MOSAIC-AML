#!/usr/bin/env python3
"""Two models x four cohorts: sensitivity and specificity, overall and per category.

  Model 1  Bulk RNA-only predictor  (BeatAML2-trained, no imputation)  -> deliverables/bulk_matrix.json
  Model 2  MOSAIC-AML multimodal    (single-cell atlas, 8 modalities)  -> scratchpad/oof_metrics.json (all-sc)
                                                                          + runs/predict_* (held-out sc)

Cohorts: BeatAML (CV), Leucegene, held-out scRNA (n=29), all scRNA (CV/OOF).
The multimodal model CANNOT run on bulk-only cohorts (BeatAML, Leucegene) - no single-cell modalities.

Every category is its own row (pipeline/categories.py says which names are one call, which
category is an umbrella over others, and which bulk category names the same lesion as a
multimodal one). Each mean says how many rows it averages, and the rows it leaves out are
listed with the reason. Each rate carries its Wilson 95% interval.

Writes, into deliverables/ (or --outdir):
  full_matrix.tsv          one row per (cell, category), plus each cell's mean
  full_matrix_paired.tsv   like-for-like: both models over the categories both report
  MOSAIC-AML_model_x_cohort_matrix.pdf/.png
Runs locally (matplotlib 3.10)."""
import argparse
import glob
import json
import math
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import categories as cat

ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
ap.add_argument("--outdir", default=os.path.join(ROOT, "deliverables"))
args = ap.parse_args()
OUTD = args.outdir; os.makedirs(OUTD, exist_ok=True)
# A cell's 95% interval wider than this is shaded in the heatmap: the value
# cannot be read to better than half the scale.
WIDE = 0.5


def entry(tp, fp, fn, tn, auroc=None):
    tp, fp, fn, tn = int(tp), int(fp), int(fn), int(tn)
    npos, nneg = tp + fn, tn + fp
    se = tp / npos if npos else float("nan")
    sp = tn / nneg if nneg else float("nan")
    (sl, sh), (pl, ph) = cat.wilson(tp, npos), cat.wilson(tn, nneg)
    return {"tp": tp, "fp": fp, "fn": fn, "tn": tn, "n_pos": npos, "n_neg": nneg,
            "sensitivity": se, "specificity": sp, "sens_lo": sl, "sens_hi": sh,
            "spec_lo": pl, "spec_hi": ph, "auroc": auroc}


def keyed(pairs, where):
    """{row: entry}, merging aliases only when they agree."""
    out = {}
    for name, e in pairs:
        k = cat.row(name)
        if k in out:
            a, b = out[k], e
            if (a["tp"], a["fp"], a["fn"], a["tn"]) != (b["tp"], b["fp"], b["fn"], b["tn"]):
                raise SystemExit("%s: %s and an alias disagree" % (where, k))
            continue
        out[k] = e
    return out


# ---- multimodal and bulk-RNA arms, all-sc (CV-OOF, same donors and folds) ----
om = json.load(open(os.path.join(ROOT, "scratchpad", "oof_metrics.json"), encoding="utf-8"))
ARMS = {a: keyed(((m, entry(v["tp"], v["fp"], v["fn"], v["tn"], v.get("auroc")))
                  for m, v in om["arms"][a]["mutations"].items()), "oof " + a)
        for a in ("multimodal", "bulkrna")}

# ---- multimodal held-out sc, from the committed reports ----
conf, counted, scores = {}, {}, {}
for f in sorted(glob.glob(os.path.join(ROOT, "runs", "predict_*", "patient_report.json"))):
    rep = json.load(open(f, encoding="utf-8"))
    for p in (rep.get("mutation_predictions") or []):
        tl = p.get("true_label")
        if tl is None:
            continue
        k = cat.row(p["mutation"])
        # an alias of a call already counted for this patient is the same call
        seen = counted.get((f, k))
        if seen is not None:
            if seen != (tl, p.get("call")):
                raise SystemExit("%s: aliases of %s disagree" % (f, k))
            continue
        counted[(f, k)] = (tl, p.get("call"))
        d = conf.setdefault(k, [0, 0, 0, 0])
        pos, cp = tl == "present", p.get("call") == "present"
        d[0 if (pos and cp) else 2 if pos else 1 if cp else 3] += 1
        if p.get("probability") is not None:
            scores.setdefault(k, ([], []))[0 if pos else 1].append(float(p["probability"]))
MM_HO = {k: entry(tp, fp, fn, tn, cat.auroc(*scores.get(k, ([], []))))
         for k, (tp, fp, fn, tn) in conf.items()}

# ---- bulk model, 4 cohorts ----
bm = json.load(open(os.path.join(ROOT, "deliverables", "bulk_matrix.json"), encoding="utf-8"))["cohorts"]


def bulk(coh):
    return keyed(((c, entry(v["tp"], v["fp"], v["fn"], v["tn"], v.get("auroc")))
                  for c, v in bm[coh]["per_category"].items()), coh)


CELLS = [("Bulk·BeatAML(CV)", "bulk", bulk("BeatAML_CV"), bm["BeatAML_CV"]["overall"]),
         ("Bulk·Leucegene", "bulk", bulk("Leucegene"), bm["Leucegene"]["overall"]),
         ("Bulk·held-out sc", "bulk", bulk("heldout_scRNA"), bm["heldout_scRNA"]["overall"]),
         ("Bulk·all sc", "bulk", bulk("all_scRNA"), bm["all_scRNA"]["overall"]),
         ("MM·held-out sc", "mm", MM_HO, None),
         ("MM·all sc(CV)", "mm", ARMS["multimodal"], None)]


def why_out(name, rows):
    """Why a row is left out of its cell's mean, or None."""
    if not rows[name]["n_pos"]:
        return "no positives among labeled patients"
    return cat.excluded(name, [k for k, v in rows.items() if v["n_pos"]])


def summary(rows):
    use = [k for k in rows if why_out(k, rows) is None]
    se = float(np.mean([rows[k]["sensitivity"] for k in use]))
    sp = float(np.nanmean([rows[k]["specificity"] for k in use]))
    au = [rows[k]["auroc"] for k in use if rows[k]["auroc"] is not None]
    return se, sp, (float(np.mean(au)) if au else None), use


SUMMARY = {}
for name, kind, rows, stored in CELLS:
    se, sp, au, use = summary(rows)
    if stored:
        # the bulk means are stored by bench_matrix.py over every category; a
        # mean over the listed rows must reproduce them
        if (len(use) != stored["n_categories"]
                or abs(se - stored["mean_sensitivity"]) > 6e-5
                or abs(sp - stored["mean_specificity"]) > 6e-5):
            raise SystemExit("%s: listed rows do not reproduce the stored mean" % name)
    SUMMARY[name] = (se, sp, au, use)

# ---- which multimodal row each bulk category matches ----
MM_NAMES = set(MM_HO) | set(ARMS["multimodal"])


def bulk_match(b):
    m, how = cat.match(b)
    if m is None and b in MM_NAMES:
        m, how = b, "same"
    return m, how


MATCHED_BY = {}
for _n, kind, rows, _s in CELLS:
    if kind == "bulk":
        for b in rows:
            m, how = bulk_match(b)
            if m:
                MATCHED_BY.setdefault(m, set()).add((b, how))


def fmt(x, nd=3):
    return "" if x is None or (isinstance(x, float) and math.isnan(x)) else "%.*f" % (nd, x)


KEYG = ["inv16", "NPM1", "FLT3", "TP53", "del7", "del5", "trisomy8", "complex", "KMT2A",
        "TET2", "DNMT3A", "ASXL1", "NRAS", "KRAS", "RUNX1", "IDH1", "IDH2", "WT1", "CEBPA",
        "KIT", "GATA2", "PTPN11"]


def order(names):
    return sorted(names, key=lambda k: (KEYG.index(cat.gene(k)) if cat.gene(k) in KEYG else len(KEYG),
                                        cat.gene(k), k))


# ---------------- per-cell TSV ----------------
COLS = ["cell", "scope", "mutation", "sensitivity", "specificity", "n_pos", "n_neg",
        "tp", "fn", "fp", "tn", "sens_lo", "sens_hi", "spec_lo", "spec_hi",
        "gene", "matches", "in_mean", "note", "auroc"]
with open(os.path.join(OUTD, "full_matrix.tsv"), "w", encoding="utf-8", newline="\n") as fh:
    fh.write("\t".join(COLS) + "\n")
    for name, kind, rows, _st in CELLS:
        se, sp, au, use = SUMMARY[name]
        fh.write("\t".join([name, "OVERALL", "(mean of %d)" % len(use), fmt(se), fmt(sp)]
                           + [""] * (len(COLS) - 6) + [fmt(au)]) + "\n")
        for k in order(rows):
            e = rows[k]
            if kind == "bulk":
                m, how = bulk_match(k)
                matches = "%s (%s)" % (m, how) if m else ""
            else:
                matches = "; ".join("%s (%s)" % bh for bh in sorted(MATCHED_BY.get(k, ())))
            out = why_out(k, rows)
            fh.write("\t".join([name, "per-mut", k, fmt(e["sensitivity"]), fmt(e["specificity"]),
                                str(e["n_pos"]), str(e["n_neg"]), str(e["tp"]), str(e["fn"]),
                                str(e["fp"]), str(e["tn"]), fmt(e["sens_lo"]), fmt(e["sens_hi"]),
                                fmt(e["spec_lo"]), fmt(e["spec_hi"]), cat.gene(k), matches,
                                "no" if out else "yes", out or "", fmt(e["auroc"])]) + "\n")

# ---------------- like-for-like ----------------
# Each comparison puts both models on the same categories. The first is the
# fairest: MOSAIC's bulk-RNA arm and the multimodal model, trained and scored
# on the same donors, folds and categories. The other two put model 1 against
# model 2 over the categories both report, on the same patients.
PAIRS = []
mm_cv = ARMS["multimodal"]
cv_use = SUMMARY["MM·all sc(CV)"][3]
PAIRS.append(("all sc (CV): bulk-RNA arm vs multimodal, same donors and folds",
              "bulk-RNA arm", "multimodal",
              [(k, k, k, "same") for k in order(cv_use) if k in ARMS["bulkrna"]],
              ARMS["bulkrna"], mm_cv))
for coh, bulk_cell, mm_cell in (("held-out sc", "Bulk·held-out sc", "MM·held-out sc"),
                                ("all sc", "Bulk·all sc", "MM·all sc(CV)")):
    brows = dict(next(c for c in CELLS if c[0] == bulk_cell)[2])
    mrows = dict(next(c for c in CELLS if c[0] == mm_cell)[2])
    shared = []
    for b in brows:
        m, how = bulk_match(b)
        if m in mrows and mrows[m]["n_pos"] and brows[b]["n_pos"]:
            shared.append((m, b, m, how))
    rank = {k: i for i, k in enumerate(order([t[0] for t in shared]))}
    shared.sort(key=lambda t: rank[t[0]])
    PAIRS.append(("%s: model 1 vs model 2, categories both report" % coh,
                  "bulk model 1", "multimodal model 2", shared, brows, mrows))

PAIRED = []


def mean_auroc(rows, names):
    """Mean AUROC over names, or None if any lacks one (held-out reports carry none)."""
    v = [rows[n]["auroc"] for n in names]
    return float(np.mean(v)) if v and all(x is not None for x in v) else None


with open(os.path.join(OUTD, "full_matrix_paired.tsv"), "w", encoding="utf-8", newline="\n") as fh:
    fh.write("comparison\tscope\tcategory\tcloseness\ta_model\ta_category\ta_sensitivity\t"
             "a_specificity\ta_auroc\ta_n_pos\tb_model\tb_category\tb_sensitivity\tb_specificity\t"
             "b_auroc\tb_n_pos\n")
    for title, an, bn, shared, arows, brows in PAIRS:
        if not shared:
            continue
        A_, B_ = [a for _k, a, _b, _h in shared], [b for _k, _a, b, _h in shared]
        a_se = float(np.mean([arows[a]["sensitivity"] for a in A_]))
        b_se = float(np.mean([brows[b]["sensitivity"] for b in B_]))
        a_sp = float(np.nanmean([arows[a]["specificity"] for a in A_]))
        b_sp = float(np.nanmean([brows[b]["specificity"] for b in B_]))
        a_au, b_au = mean_auroc(arows, A_), mean_auroc(brows, B_)
        PAIRED.append((title, an, bn, len(shared), a_se, b_se, a_sp, b_sp, a_au, b_au))
        fh.write("\t".join([title, "OVERALL", "(mean of %d)" % len(shared), "", an, "",
                            fmt(a_se), fmt(a_sp), fmt(a_au), "", bn, "", fmt(b_se), fmt(b_sp),
                            fmt(b_au), ""]) + "\n")
        for k, a, b, how in shared:
            ea, eb = arows[a], brows[b]
            fh.write("\t".join([title, "per-mut", k, how, an, a, fmt(ea["sensitivity"]),
                                fmt(ea["specificity"]), fmt(ea["auroc"]), str(ea["n_pos"]), bn, b,
                                fmt(eb["sensitivity"]), fmt(eb["specificity"]), fmt(eb["auroc"]),
                                str(eb["n_pos"])]) + "\n")

# ---------------- figure ----------------
C_SENS, C_SPEC = "#2a78d6", "#008300"         # validated categorical slots 1 & 2
INK, MUTED, GRID, SURF, NA = "#1a1a19", "#6b6b63", "#d9d9d4", "#ffffff", "#eceae2"
names = [c[0] for c in CELLS]

# heatmap rows: a bulk category that names a multimodal row's lesion shares its row
ROWS = {}                                      # row label -> {cell: (entry, category)}
for name, kind, rows, _st in CELLS:
    for k, e in rows.items():
        if not e["n_pos"]:
            continue
        if kind == "bulk":
            m, how = bulk_match(k)
            label = k if not m or m == k else "%s %s %s" % (m, "=" if how == "same" else "≈", k)
        else:
            bs = sorted(MATCHED_BY.get(k, ()))
            label = k if not bs or bs[0][0] == k else "%s %s %s" % (k, "=" if bs[0][1] == "same" else "≈", bs[0][0])
        out = why_out(k, rows)
        if out:
            label += "  (not in mean)"
        ROWS.setdefault(label, {})[name] = (e, k)
row_key = {lab: next(iter(v.values()))[1] for lab, v in ROWS.items()}
rank = {k: i for i, k in enumerate(order(set(row_key.values())))}
labels = sorted(ROWS, key=lambda lab: (rank[row_key[lab]], lab))

nrow = len(labels)
fig = plt.figure(figsize=(14.5, 7.6 + 0.24 * nrow), facecolor=SURF)
gs = fig.add_gridspec(3, 1, height_ratios=[3.4, 2.6, 0.24 * nrow], hspace=0.22)

# -- panel A: each cell's own mean --
axA = fig.add_subplot(gs[0]); axA.set_facecolor(SURF)
x = np.arange(len(names)); w = 0.38
se = [SUMMARY[n][0] for n in names]; sp = [SUMMARY[n][1] for n in names]
au = [SUMMARY[n][2] for n in names]; nn = [len(SUMMARY[n][3]) for n in names]
axA.bar(x - w/2, se, w, color=C_SENS, edgecolor=SURF, lw=1.5, label="sensitivity")
axA.bar(x + w/2, sp, w, color=C_SPEC, edgecolor=SURF, lw=1.5, label="specificity")
for i in range(len(names)):
    axA.text(x[i]-w/2, se[i]+0.02, "%.2f" % se[i], ha="center", va="bottom", fontsize=8, color=C_SENS, fontweight="bold")
    axA.text(x[i]+w/2, sp[i]+0.02, "%.2f" % sp[i], ha="center", va="bottom", fontsize=8, color=C_SPEC, fontweight="bold")
    axA.text(x[i], -0.12, "mean of %d" % nn[i] + ("\nAUROC %.2f" % au[i] if au[i] is not None else ""),
             ha="center", va="top", fontsize=7.5, color=MUTED)
for xi in (len(names), len(names) + 1):
    axA.axvspan(xi-0.5, xi+0.5, color=NA, zorder=0)
    axA.text(xi, 0.5, "N/A\nbulk cohort:\nno single-cell\nmodalities", ha="center", va="center",
             fontsize=8, color=MUTED, style="italic")
axA.set_xlim(-0.6, len(names)+1.6); axA.set_ylim(0, 1.06)
axA.set_xticks(list(x)+[len(names), len(names)+1])
axA.set_xticklabels(names+["MM·BeatAML", "MM·Leucegene"], fontsize=8.7)
axA.set_ylabel("score", fontsize=10, color=INK); axA.set_yticks([0, .25, .5, .75, 1.0])
axA.spines[["top", "right"]].set_visible(False); axA.spines[["left", "bottom"]].set_color(GRID)
axA.tick_params(colors=MUTED, labelsize=8.5); axA.yaxis.grid(True, color=GRID, lw=.7); axA.set_axisbelow(True)
axA.axvline(3.5, color=GRID, lw=1, ls=(0, (4, 3)))
axA.text(1.5, 1.10, "Bulk RNA-only model ①  (BeatAML-trained)", ha="center", fontsize=9.5, color=INK, fontweight="bold")
axA.text(5.4, 1.10, "MOSAIC-AML multimodal ②  (single-cell)", ha="center", fontsize=9.5, color=INK, fontweight="bold")
axA.set_title("A  Each cell's mean over its own categories (the categories differ between cells)",
              fontsize=10.5, color=INK, loc="left", pad=22)

# -- panel B: like-for-like --
axB = fig.add_subplot(gs[1]); axB.set_facecolor(SURF)
for i, (title, an, bn, k, a_se, b_se, a_sp, b_sp, a_au, b_au) in enumerate(PAIRED):
    base = i * 5.0
    for j, (v, col, lab) in enumerate([(a_se, C_SENS, an), (b_se, C_SENS, bn), (a_sp, C_SPEC, an), (b_sp, C_SPEC, bn)]):
        xx = base + j + (0.4 if j >= 2 else 0)
        axB.bar(xx, v, 0.8, color=col, alpha=(0.55 if j % 2 == 0 else 1.0), edgecolor=SURF)
        axB.text(xx, v + 0.02, "%.2f" % v, ha="center", va="bottom", fontsize=7.5, color=INK)
        axB.text(xx, -0.04, lab.replace(" model ", "\nmodel ").replace(" arm", "\narm"), ha="center", va="top", fontsize=6.6, color=MUTED)
    au_txt = ("; mean AUROC %.2f vs %.2f" % (a_au, b_au)) if a_au is not None and b_au is not None else ""
    axB.text(base + 1.7, 1.12, "%s\n(mean of %d%s)" % (title, k, au_txt), ha="center", va="bottom", fontsize=7.6, color=INK)
axB.set_xlim(-0.8, 5.0 * len(PAIRED) - 0.8); axB.set_ylim(0, 1.06)
axB.set_xticks([]); axB.set_yticks([0, .5, 1.0])
axB.spines[["top", "right", "bottom"]].set_visible(False); axB.spines["left"].set_color(GRID)
axB.tick_params(colors=MUTED, labelsize=8); axB.yaxis.grid(True, color=GRID, lw=.7); axB.set_axisbelow(True)
axB.set_title("B  Like-for-like: both models over the same categories  (blue sensitivity, green specificity; pale = first model)",
              fontsize=10.5, color=INK, loc="left", pad=30)

# -- panel C: per-category sensitivity --
axC = fig.add_subplot(gs[2]); axC.set_facecolor(SURF)
M = np.full((nrow, len(CELLS)), np.nan)
for i, lab in enumerate(labels):
    for j, name in enumerate(names):
        if name in ROWS[lab]:
            M[i, j] = ROWS[lab][name][0]["sensitivity"]
im = axC.imshow(M, aspect="auto", cmap="YlGnBu", vmin=0, vmax=1)
axC.set_xticks(range(len(CELLS))); axC.set_xticklabels(names, fontsize=8.5, rotation=20, ha="right")
axC.set_yticks(range(nrow)); axC.set_yticklabels(labels, fontsize=7.6)
for i, lab in enumerate(labels):
    for j, name in enumerate(names):
        if name not in ROWS[lab]:
            axC.add_patch(plt.Rectangle((j-0.5, i-0.5), 1, 1, color=NA, zorder=1))
            continue
        e = ROWS[lab][name][0]
        wide = e["sens_hi"] - e["sens_lo"] > WIDE
        if wide:
            # faded: the value is too uncertain to read as a color
            axC.add_patch(plt.Rectangle((j-0.5, i-0.5), 1, 1, color=SURF, alpha=0.62, lw=0, zorder=2))
        axC.text(j, i, "%.2f  %d/%d" % (e["sensitivity"], e["tp"], e["n_pos"]), ha="center", va="center",
                 fontsize=6.4, color=(MUTED if wide else "white" if M[i, j] > 0.55 else INK),
                 style=("italic" if wide else "normal"), zorder=3)
axC.xaxis.tick_top()
axC.set_xticklabels(names, fontsize=8.5, rotation=0, ha="center")
axC.set_title("C  Per-category sensitivity, with TP/positives.  Faded: 95%% interval wider than %.1f.  "
              "\"=\" same lesion, \"≈\" nearly the same.  Grey: not reported." % WIDE,
              fontsize=10.0, color=INK, loc="left", pad=24)
cb = fig.colorbar(im, ax=axC, fraction=0.02, pad=0.01); cb.set_label("sensitivity", fontsize=8, color=MUTED)
cb.ax.tick_params(labelsize=7, colors=MUTED)

handles = [Patch(facecolor=C_SENS, label="sensitivity"), Patch(facecolor=C_SPEC, label="specificity")]
fig.legend(handles=handles, loc="upper center", ncol=2, frameon=False, fontsize=9.5,
           bbox_to_anchor=(0.5, 0.985), labelcolor=INK)
fig.suptitle("MOSAIC-AML - mutation-calling performance by model and cohort", fontsize=15.5,
             color=INK, fontweight="bold", y=0.998)
fig.text(0.5, 0.004,
         "Bulk model ① is BeatAML2-trained (bulk RNA, no imputation); multimodal ② is single-cell-trained (8 modalities). "
         "Both operate at the F1-max threshold; in all-sc CV that threshold is chosen on the vector it is scored on, so those "
         "absolute rates are optimistic and only contrasts between arms are fair. Rows a mean leaves out: an umbrella category "
         "listed beside its subtypes, a second model head for one lesion, a category with no positives.",
         ha="center", fontsize=7.0, color=MUTED, style="italic", wrap=True)
fig.subplots_adjust(left=0.24, right=0.97, top=0.93, bottom=0.025)
for ext, dpi in [("pdf", None), ("png", 130)]:
    fig.savefig(os.path.join(OUTD, "MOSAIC-AML_model_x_cohort_matrix.%s" % ext),
                **({"dpi": dpi} if dpi else {}), facecolor=SURF)
print("wrote %s: full_matrix.tsv, full_matrix_paired.tsv, MOSAIC-AML_model_x_cohort_matrix.pdf/.png" % OUTD)
for title, an, bn, k, a_se, b_se, a_sp, b_sp, a_au, b_au in PAIRED:
    print("  %-70s k=%2d  sens %.3f vs %.3f  spec %.3f vs %.3f  auroc %s vs %s"
          % (title, k, a_se, b_se, a_sp, b_sp, fmt(a_au), fmt(b_au)))
print("MATRIX OK")
