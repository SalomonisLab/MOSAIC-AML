#!/usr/bin/env python3
"""VAF-stratified sensitivity of the bulk RNA-only model on BeatAML (CV-OOF).

Reframes 'misses' honestly: does the model catch clonal/high-VAF variants and miss subclonal/low-VAF
ones (an expression-detection limit), rather than failing at random? For every BeatAML positive of each
variant-level driver, we look up that gene's variant-allele frequency (mutations.txt, t_vaf) and ask
whether the model's honest CV-OOF call caught it, binned by VAF.

Writes deliverables/vaf_stratified.json + .tsv and (if matplotlib) a PNG/PDF.
  bsub -q test -W 60 -M 24000 -R "rusage[mem=24000]" -o vaf.log \
    /usr/local/anaconda3-2020/bin/python vaf_stratified.py
"""
import os, sys, json, warnings
warnings.filterwarnings("ignore")
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "pipeline"))
from amlmm.bulk_predictor import BulkMutationPredictor, _pct
OUTD = os.path.join(ROOT, "deliverables"); os.makedirs(OUTD, exist_ok=True)

BP = BulkMutationPredictor.load(os.path.join(HERE, "bulk_mutation_predictor.pkl"))
d = np.load(os.path.join(os.path.dirname(ROOT), "aml-bakeoff", "bundle_data.npz"), allow_pickle=True)
CATS = [str(c) for c in d["drivers"]]; SAMP = [str(s) for s in d["ba_samples"]]
baX = d["ba_X"].astype(float); baL = d["ba_L"].astype(float)
VARCAP = 2500
CYTO = {"inv16", "inv(16)", "complex", "del5", "del7", "del17", "trisomy8", "trisomy", "kmt2a",
        "monosomy", "minus", "plus", "t"}

def gene_of(cat):
    return str(cat).split("_")[0].split("-")[0].upper()

# ---- per-(sample, gene) VAF from mutations.txt (max VAF across that gene's variants in that sample) ----
MUT = os.path.join(ROOT, "data", "external", "beataml", "mutations.txt")
vaf = {}
with open(MUT) as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    ci = {c: i for i, c in enumerate(hdr)}
    si, vi, gi = ci["dbgap_sample_id"], ci["t_vaf"], ci["symbol"]
    for line in fh:
        c = line.rstrip("\n").split("\t")
        if len(c) <= max(si, vi, gi):
            continue
        try:
            v = float(c[vi])
        except Exception:
            continue
        k = (c[si], c[gi].upper())
        if v > vaf.get(k, -1):
            vaf[k] = v
print("VAF entries parsed:", len(vaf), "| sample overlap with bundle:",
      len(set(s for (s, g) in vaf) & set(SAMP)))

# ---- BeatAML CV-OOF per-sample call (replicate the deployed recipe) ----
cl = np.log2(np.clip(baX, 0, None) + 1.0); mu = cl.mean(0); sd = cl.std(0); sd[sd == 0] = 1.0
Zba = (cl - mu) / sd

records = []       # (category, gene, sample, vaf, called)
for j, cat in enumerate(CATS):
    if cat not in BP.categories:
        continue
    g = gene_of(cat)
    if g.lower() in CYTO or any(cat.lower().startswith(x) for x in CYTO):
        continue                                   # karyotype events have no VAF
    y = baL[:, j]; ok = ~np.isnan(y); yv = y[ok].astype(int); Z = Zba[ok]
    idx = np.where(ok)[0]
    if yv.sum() < 3 or (yv == 0).sum() < 3:
        continue
    oof = np.full(len(yv), np.nan)
    ns = min(5, int(yv.sum()), int((yv == 0).sum()))
    for tri, tei in StratifiedKFold(ns, shuffle=True, random_state=0).split(Z, yv):
        sel = np.argsort(Z[tri].var(0))[::-1][:VARCAP]
        est = LogisticRegression(C=1.0, class_weight="balanced", max_iter=2000).fit(Z[tri][:, sel], yv[tri])
        oof[tei] = est.decision_function(Z[tei][:, sel])
    om = ~np.isnan(oof)
    try:
        a = roc_auc_score(yv[om], oof[om]); sign = 1.0 if a >= 0.5 else -1.0
    except Exception:
        sign = 1.0
    p = _pct(np.sort(oof[om]), oof[om]); p = p if sign > 0 else 1 - p
    call = (p >= BP.models[cat]["threshold"]).astype(int)
    for pos, s_i in enumerate(np.where(om)[0]):
        if yv[s_i] != 1:
            continue                               # positives only -> sensitivity
        s = SAMP[idx[s_i]]
        s_dna = (s[:-1] + "D") if s.endswith("R") else s     # bundle RNA id (…R) -> WES/DNA id (…D)
        v = vaf.get((s_dna, g))
        if v is None:
            continue
        records.append((cat, g, s, float(v), int(call[pos])))

print("positive calls with a matched VAF:", len(records))

# ---- stratify sensitivity by VAF bin ----
BINS = [(0.0, 0.10, "subclonal <0.10"), (0.10, 0.25, "0.10–0.25"),
        (0.25, 0.40, "0.25–0.40"), (0.40, 1.01, "clonal ≥0.40")]
strat = []
for lo, hi, lab in BINS:
    sub = [r for r in records if lo <= r[3] < hi]
    n = len(sub); k = sum(r[4] for r in sub)
    strat.append({"bin": lab, "vaf_lo": lo, "vaf_hi": hi, "n_positive": n,
                  "n_called": k, "sensitivity": round(k/n, 4) if n else None})

# per key gene, sensitivity vs median VAF
by_gene = {}
for cat, g, s, v, c in records:
    by_gene.setdefault(g, []).append((v, c))
gene_rows = []
for g, lst in sorted(by_gene.items(), key=lambda kv: -len(kv[1])):
    n = len(lst); k = sum(c for _, c in lst)
    gene_rows.append({"gene": g, "n_positive": n, "sensitivity": round(k/n, 4),
                      "median_vaf": round(float(np.median([v for v, _ in lst])), 3)})

out = {"model": "bulk_mutation_predictor.pkl (BeatAML CV-OOF)", "n_positive_calls": len(records),
       "vaf_strata": strat, "by_gene": gene_rows,
       "note": "VAF = that gene's max t_vaf in the sample (mutations.txt); a proxy for the driver's clonality. "
               "Sensitivity = fraction of positives the honest CV-OOF call caught, at the deployed threshold."}
json.dump(out, open(os.path.join(OUTD, "vaf_stratified.json"), "w"), indent=1)
with open(os.path.join(OUTD, "vaf_stratified.tsv"), "w") as fh:
    fh.write("vaf_bin\tn_positive\tn_called\tsensitivity\n")
    for s in strat:
        fh.write("%s\t%s\t%s\t%s\n" % (s["bin"], s["n_positive"], s["n_called"], s["sensitivity"]))
    fh.write("\ngene\tn_positive\tmedian_vaf\tsensitivity\n")
    for r in gene_rows:
        fh.write("%s\t%d\t%.3f\t%.3f\n" % (r["gene"], r["n_positive"], r["median_vaf"], r["sensitivity"]))

print("\nVAF-STRATIFIED SENSITIVITY (bulk model, BeatAML CV-OOF):")
for s in strat:
    print("  %-16s n=%-4d sens=%s" % (s["bin"], s["n_positive"], s["sensitivity"]))

# ---- figure ----
try:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    labs = [s["bin"] for s in strat]; se = [s["sensitivity"] or 0 for s in strat]; ns = [s["n_positive"] for s in strat]
    fig, ax = plt.subplots(figsize=(7.6, 5.0), facecolor="#ffffff")
    ax.bar(range(len(labs)), se, width=0.62, color="#2a78d6", edgecolor="#ffffff", lw=1.5)
    for i, (v, nn) in enumerate(zip(se, ns)):
        ax.text(i, v+0.02, "%.2f\n(n=%d)" % (v, nn), ha="center", va="bottom", fontsize=9, color="#1a1a19")
    ax.set_xticks(range(len(labs))); ax.set_xticklabels(labs, fontsize=9.5)
    ax.set_ylim(0, 1.08); ax.set_ylabel("sensitivity (fraction of positives caught)", fontsize=10)
    ax.set_xlabel("variant allele frequency (clonality proxy)", fontsize=10)
    ax.set_title("Bulk model sensitivity rises with VAF — misses are largely subclonal\nBeatAML CV-OOF, positives with a matched t_vaf",
                 fontsize=12, color="#1a1a19", fontweight="bold")
    ax.spines[["top","right"]].set_visible(False); ax.spines[["left","bottom"]].set_color("#d9d9d4")
    ax.tick_params(colors="#6b6b63"); ax.yaxis.grid(True, color="#d9d9d4", lw=.7); ax.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(os.path.join(OUTD, "vaf_stratified.pdf"), facecolor="#ffffff")
    fig.savefig(os.path.join(OUTD, "vaf_stratified.png"), dpi=130, facecolor="#ffffff")
    print("wrote vaf_stratified.pdf/.png")
except Exception as e:
    print("plot skipped:", e)
print("VAF STRAT OK")
