#!/usr/bin/env python3
"""Validate the MOSAIC-AML mutation predictor against the Trumpp/Waclawiczek cohort known mutations.

Trumpp samples are NEW EXTERNAL, so only COMPOSITION is available (from the cellHarmony alignment) —
the imputed modalities (ADT/Lipid/Metabolite/GRN) that carry most of the predictor's signal were never
computed for them. So this is the HONEST degraded composition-mode test, NOT the ~0.86 cohort number.

Builds each sample's 89-state composition from cellHarmony_lite_assignments.txt, runs the predictor
(composition only), and compares to the Table S4 known drivers. Reports per-mutation AUC across the 16
samples + writes board reports.  Run on an LSF compute node.
"""
import os, sys, json, pickle, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
from amlmm.predictor import MutationPredictor
try:
    import control_gate as CG
except Exception:
    CG = None

HERE = os.path.dirname(os.path.abspath(__file__))
ASSIGN = "/data/salomonis2/LabFiles/Frank-Li/scTriangulate/Hs_AML_UDON/output/cellHarmony_lite_assignments.txt"
RUNS_ROOT = "/data/salomonis-archive/LabFiles/Nicholas/AML-multimodal/runs"
def log(m): print(m, flush=True)

# Table S4 ground truth -> predictor flags (NF1/ASXL1 noted; del5/kmt2a from karyotype). del5=del(5), +8=trisomy8.
KNOWN = {
 "P9_Diagnosis.Heidelberg":  {"complex"},
 "P9_Refractory.Heidelberg": {"del5"},
 "P10_Diagnosis.Heidelberg": {"TP53","NRAS","complex"},
 "P10_Refractory.Heidelberg":{"TP53","NRAS","complex"},
 "P11_Diagnosis.Heidelberg": {"ASXL1","complex"},
 "P11_Refractory.Heidelberg":{"ASXL1","complex"},
 "P13_Diagnosis.Heidelberg": {"NPM1","FLT3","DNMT3A"},
 "P13_Refractory.Heidelberg":{"NPM1","FLT3"},
 "P14_Diagnosis.Heidelberg": {"ASXL1","FLT3","NPM1","TET2","trisomy8"},
 "P14_Refractory.Heidelberg":{"ASXL1","FLT3","NPM1","TET2","trisomy8"},
 "P16_Diagnosis.Heidelberg": {"kmt2a"},
 "P16_Refractory.Heidelberg":{"TP53","kmt2a"},
 "P18_Diagnosis.Heidelberg": {"FLT3","IDH2","NPM1"},
 "P18_Refractory.Heidelberg":{"FLT3","IDH2","NPM1"},
 "P19_Diagnosis.Heidelberg": {"DNMT3A","TET2","TP53","complex"},
 "P19_Refractory.Heidelberg":{"DNMT3A","TET2","TP53","complex"},
}

P = MutationPredictor.load(os.path.join(HERE, "mutation_predictor.pkl"))
log("predictor: %d mutations, modalities %s" % (len(P.mutations), P.modalities))
comp_features = None
for (mut, mod), mm in P.models.items():
    if mod == "Composition":
        comp_features = mm.features; break
log("composition features: %d" % (len(comp_features) if comp_features else 0))

# ---- per-sample composition from cellHarmony assignments ----
counts = {}
with open(ASSIGN, encoding="utf-8", errors="replace") as fh:
    header = fh.readline().rstrip("\n").split("\t")
    bc_i, st_i = header.index("CellBarcode"), header.index("Hs-BM-titrated-reference-centroid")
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) <= max(bc_i, st_i):
            continue
        samp = p[bc_i].split(".", 1)[1] if "." in p[bc_i] else None
        if samp is None:
            continue
        counts.setdefault(samp, {}).setdefault(p[st_i], 0)
        counts[samp][p[st_i]] += 1
comp = {}
for s, d in counts.items():
    tot = sum(d.values()) or 1
    comp[s] = pd.Series({k: v / tot for k, v in d.items()})
log("samples with composition: %d (%s)" % (len(comp), ", ".join(sorted(comp)[:3]) + " ..."))

# ---- predict (composition only) + collect for AUC ----
gate = CG.load_gate() if CG else None
scope = [m for m in P.mutations if any((m in KNOWN[s]) for s in KNOWN)]    # mutations that appear in the truth
y_by_mut = {m: [] for m in scope}; s_by_mut = {m: [] for m in scope}
rows = []
for s in sorted(KNOWN):
    if s not in comp:
        log("  no composition for %s" % s); continue
    sample = {"Composition": comp[s]}
    preds = []
    for m in P.mutations:
        if (m, "Composition") not in P.models:
            continue
        pr = P.predict_one(m, sample); pr["mutation"] = m
        truth = "present" if m in KNOWN[s] else "absent"
        pr["true_label"] = truth
        preds.append(pr)
        if m in scope and pr["probability"] is not None:
            y_by_mut[m].append(1 if m in KNOWN[s] else 0); s_by_mut[m].append(pr["probability"])
    preds.sort(key=lambda p: -(p["probability"] or 0))
    spec = CG.score_gate(gate, comp[s]) if gate is not None else None
    rep = {"mode": "mutation_panel", "sample_key": "Trumpp::" + s, "dataset": "Trumpp/Waclawiczek",
           "specimen_class": (spec["call"] if spec else None), "control_gate": spec,
           "mutation_predictions": preds, "predictor": P.summary(), "validation": True,
           "modalities_available": ["Composition"], "known_drivers": sorted(KNOWN[s]),
           "note": "external Trumpp sample — COMPOSITION-ONLY (imputed modalities not available); predicted vs Table S4"}
    d = os.path.join(RUNS_ROOT, "trumpp_" + "".join(ch if ch.isalnum() else "_" for ch in s))
    os.makedirs(d, exist_ok=True)
    json.dump(rep, open(os.path.join(d, "patient_report.json"), "w"), default=str, indent=1)
    # per-sample: did the top predicted present-calls match known?
    called = [p["mutation"] for p in preds if p["call"] == "present"]
    rows.append((s, sorted(KNOWN[s] & set(P.mutations)), called))

log("\n=== Trumpp per-mutation AUC (composition-only, across %d samples) ===" % len(KNOWN))
log("%-10s %5s %5s %7s" % ("mutation", "n+", "n", "AUC"))
aucs = []
for m in scope:
    y, sc = np.array(y_by_mut[m]), np.array(s_by_mut[m])
    if len(set(y)) == 2 and len(y) >= 6:
        a = roc_auc_score(y, sc); aucs.append(a)
        log("%-10s %5d %5d %7.3f" % (m, int(y.sum()), len(y), a))
    else:
        log("%-10s %5d %5d   n/a (need both classes)" % (m, int(y.sum()), len(y)))
log("MEAN AUC = %.3f over %d evaluable mutations" % (np.mean(aucs) if aucs else float("nan"), len(aucs)))

log("\n=== per-sample: known drivers (in panel) vs called-present ===")
for s, known, called in rows:
    log("%-26s known=%-28s called=%s" % (s, ",".join(known) or "—", ",".join(called) or "—"))
log("\nwrote %d Trumpp board reports -> %s/trumpp_*" % (len(rows), RUNS_ROOT))
log("VALIDATE TRUMPP OK")
