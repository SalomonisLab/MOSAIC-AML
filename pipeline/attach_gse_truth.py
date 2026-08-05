#!/usr/bin/env python3
"""Attach GSE281087 ground truth to the ingested reports so they can be evaluated.

`ingest_patient.py` is the new-patient path and does not know any genotype, so its reports carry
`true_label: null`. GSE281087 ships a Mutation_Matrix AND a Panel_Coverage sheet, so we can label
honestly: a gene the panel ASSAYED gets present/absent; a gene it NEVER assayed is left unlabelled
(null) and therefore excluded from the metrics — a 0 there means "not tested", not wild-type.

  /usr/local/anaconda3-2020/bin/python attach_gse_truth.py    (runs locally too)
"""
import os, csv, json, glob
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
P = os.environ.get("GSE_DIR", "/data/salomonis-archive/FASTQs/PublicDatasets/DATASET_8_GSE281087")
MMP = os.path.join(P, "metadata_harmonized", "GSE281087_Mutation_Matrix.csv")
PCP = os.path.join(P, "metadata_harmonized", "GSE281087_Panel_Coverage.csv")
if not os.path.exists(MMP):                                     # local fallback copy
    alt = os.path.join(ROOT, "labels", "GSE281087_Mutation_Matrix.csv")
    if os.path.exists(alt): MMP, PCP = alt, os.path.join(ROOT, "labels", "GSE281087_Panel_Coverage.csv")

rows = list(csv.reader(open(MMP)))
hdr = rows[0]; gcols = hdr[2:]
truth = {}
for r in rows[1:]:
    truth[r[1].strip()] = {gcols[i].upper(): r[2 + i].strip() for i in range(len(gcols))}
assayed = {}
for r in csv.DictReader(open(PCP)):
    assayed[r["template_column"].upper()] = (r["assayed_in_GSE281087"].strip().lower() == "yes")
print("truth samples %d | genes assayed by the panel %d/%d" % (len(truth), sum(assayed.values()), len(assayed)))

def gene_of(c):
    cl = str(c).lower()
    if "inv16" in cl or "inv(16)" in cl or "cbfb" in cl: return "INV16"
    if "kmt2a" in cl: return "KMT2A"
    u = str(c).upper()
    if u.startswith("FLT3_ITD") or "ITD" in u: return "FLT3-ITD"
    if u.startswith("FLT3_TKD") or "TKD" in u: return "FLT3-TKD"
    return str(c).split("_")[0].split("-")[0].upper()

n_rep = n_lab = n_skip = 0
for f in sorted(glob.glob(os.path.join(ROOT, "runs", "gse_*", "patient_report.json"))):
    rep = json.load(open(f))
    key = str(rep.get("sample_key", ""))
    name = key.split("::")[-1]
    t = truth.get(name)
    if t is None:                                                # slug/underscore variants
        cand = [s for s in truth if s.replace("-", "_") == name.replace("-", "_")]
        t = truth.get(cand[0]) if cand else None
    if t is None:
        n_skip += 1; continue
    labelled = 0
    for p in (rep.get("mutation_predictions") or []):
        g = gene_of(p["mutation"])
        if not assayed.get(g, False):                            # NEVER assayed -> leave unlabelled
            p["true_label"] = None; continue
        v = t.get(g)
        if v in ("0", "1"):
            p["true_label"] = "present" if v == "1" else "absent"; labelled += 1
        else:
            p["true_label"] = None
    rep["validation"] = True
    rep["truth_source"] = "GSE281087_Mutation_Matrix (panel-honest: never-assayed genes left unlabelled)"
    json.dump(rep, open(f, "w"), default=str, indent=1)
    n_rep += 1; n_lab += labelled
print("labelled %d reports, %d gene calls total (%d reports had no matching truth row)" % (n_rep, n_lab, n_skip))
print("ATTACH GSE TRUTH OK")
