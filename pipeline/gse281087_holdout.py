#!/usr/bin/env python3
"""New external single-cell held-out cohort: GSE281087 (~15 scRNA samples with mutation calls).

Each sample's 10x .h5 → bulk-equivalent (sum cells → CP10k) → scored by the BULK RNA-only caller
(ref='sc'), gene-level, compared to GSE281087_Mutation_Matrix. Crucially uses GSE281087_Panel_Coverage
to compute HONEST specificity: a gene the panel never assayed is EXCLUDED from the negatives (a `0`
there is "not tested", not wild-type) — realizing the genes-assayed correction on a real cohort.

Writes deliverables/gse281087_holdout.json + .tsv.
  bsub -q test -W 30 -M 12000 -R "rusage[mem=12000]" -o gse.log \
    /usr/local/anaconda3-2020/bin/python gse281087_holdout.py
"""
import os, sys, glob, csv, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, h5py
from amlmm.bulk_predictor import BulkMutationPredictor

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
OUTD = os.path.join(ROOT, "deliverables"); os.makedirs(OUTD, exist_ok=True)
P = "/data/salomonis-archive/FASTQs/PublicDatasets/DATASET_8_GSE281087"
BP = BulkMutationPredictor.load(os.path.join(HERE, "bulk_mutation_predictor.pkl"))

def gene_of(cat): return str(cat).split("_")[0].split("-")[0].upper()

# ---- truth + panel coverage ----
mm_rows = list(csv.reader(open(P+"/metadata_harmonized/GSE281087_Mutation_Matrix.csv")))
hdr = mm_rows[0]; gcols = hdr[2:]
truth = {}                                    # sample -> {GENE: 0/1}
for r in mm_rows[1:]:
    truth[r[1]] = {gcols[i].upper(): int(r[2+i]) for i in range(len(gcols)) if r[2+i] in ("0","1")}
assayed = {}                                  # GENE -> assayed? (panel coverage)
for r in csv.DictReader(open(P+"/metadata_harmonized/GSE281087_Panel_Coverage.csv")):
    assayed[r["template_column"].upper()] = (r["assayed_in_GSE281087"].strip().lower() == "yes")
print("truth samples: %d | panel genes assayed: %d/%d" %
      (len(truth), sum(assayed.values()), len(assayed)))

# ---- read a 10x .h5 -> per-gene CP10k bulk-equivalent ----
def bulk_equiv(h5path):
    with h5py.File(h5path, "r") as h:
        m = h["matrix"]; data = m["data"][:]; idx = m["indices"][:]
        names = [x.decode() if isinstance(x, bytes) else str(x) for x in m["features"]["name"][:]]
        ng = len(names)
        gsum = np.bincount(idx, weights=data, minlength=ng).astype(float)   # indices = gene row -> per-gene total counts
    tot = gsum.sum() or 1.0
    import pandas as pd
    return pd.Series(gsum / tot * 1e4, index=names)                        # CP10k linear

# ---- score each sample, gene-level calls ----
import pandas as pd
rows, per_sample = [], {}
for f in sorted(glob.glob(P+"/counts_h5/*.h5")):
    name = os.path.basename(f).split("_counts.h5")[0].split("_", 1)[1]     # GSM..._JH4189_single -> JH4189_single
    if name not in truth:
        alt = [s for s in truth if s.split("_")[0] == name.split("_")[0]]  # match on donor prefix if suffix differs
        if not alt:
            print("  no truth for", name); continue
        name = alt[0]
    ser = bulk_equiv(f)
    z = BP._z(BP._clog(BP._align(ser)), "sc")
    call_gene = {}
    for c in BP.categories:
        g = gene_of(c)
        pr = BP.predict_one(c, z, "sc")
        if pr["call"] == "present":
            call_gene[g] = 1
        else:
            call_gene.setdefault(g, 0)
    per_sample[name] = call_gene

# ---- confusion, panel-honest vs naive ----
def score(restrict_assayed):
    tp=fp=fn=tn=0; pg={}
    for s, calls in per_sample.items():
        t = truth.get(s, {})
        for g in set(list(calls) + list(t)):
            G = g.upper()
            if G not in truth[s] and G not in calls: continue
            if restrict_assayed and not assayed.get(G, False): continue     # skip never-assayed genes
            if G not in BP_GENES: continue                                  # only genes the caller can predict
            y = t.get(G, 0); c = calls.get(G, 0)
            d = pg.setdefault(G, [0,0,0,0])
            if y and c: tp+=1; d[0]+=1
            elif y and not c: fn+=1; d[2]+=1
            elif (not y) and c: fp+=1; d[1]+=1
            else: tn+=1; d[3]+=1
    se = tp/(tp+fn) if tp+fn else None; sp = tn/(tn+fp) if tn+fp else None
    return dict(tp=tp,fp=fp,fn=fn,tn=tn,sensitivity=round(se,4) if se is not None else None,
                specificity=round(sp,4) if sp is not None else None), pg

BP_GENES = set(gene_of(c) for c in BP.categories)
honest, pg_h = score(True); naive, pg_n = score(False)
out = {"cohort":"GSE281087 (single-cell external held-out)","n_samples_scored":len(per_sample),
       "n_truth_samples":len(truth),"caller":"bulk_mutation_predictor.pkl (bulk-equivalent, ref=sc)",
       "panel_honest":honest,"naive_all_genes":naive,
       "note":"panel_honest excludes never-assayed genes from the negatives (0 != wild-type there); "
              "naive counts them as negatives and so UNDER-states specificity."}
json.dump(out, open(os.path.join(OUTD,"gse281087_holdout.json"),"w"), indent=1)
with open(os.path.join(OUTD,"gse281087_holdout.tsv"),"w") as fh:
    fh.write("scope\tTP\tFP\tFN\tTN\tsensitivity\tspecificity\n")
    for lab,m in [("panel-honest",honest),("naive-all-genes",naive)]:
        fh.write("%s\t%d\t%d\t%d\t%d\t%s\t%s\n"%(lab,m["tp"],m["fp"],m["fn"],m["tn"],m["sensitivity"],m["specificity"]))
print("samples scored:", len(per_sample))
print("PANEL-HONEST : sens=%s spec=%s (TP=%d FP=%d FN=%d TN=%d)"%(honest["sensitivity"],honest["specificity"],honest["tp"],honest["fp"],honest["fn"],honest["tn"]))
print("naive(all)   : sens=%s spec=%s (FP=%d)"%(naive["sensitivity"],naive["specificity"],naive["fp"]))
print("GSE281087 HOLDOUT OK")
