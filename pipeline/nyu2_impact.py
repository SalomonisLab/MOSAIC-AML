#!/usr/bin/env python3
"""Impact of the new NYU-2 metadata (AML_harmonized_metadata_v2_NYU2.xlsx).

NYU-2's 28 samples previously had NO mutation calls in mutation_matrix_explicit_v2.tsv — every gene read
as 0, i.e. they were silently scored as WILD-TYPE for everything. The new file supplies real calls. This
is exactly the "absent can mean not-assayed" label defect documented in truth_provenance.md, now fixable.

Reports, per mutation: current n_pos, NYU-2 positives added (restricted to samples that actually have
atlas expression), new n_pos, and which mutations cross the >=8-positive trainability threshold.

  bsub -q test -W 20 -M 16000 -R "rusage[mem=16000]" -o nyu2.log \
    /usr/local/anaconda3-2020/bin/python nyu2_impact.py <path-to-xlsx>
"""
import os, sys, csv, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pandas as pd
from amlmm.context import build_context, Config

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
XLSX = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "labels", "AML_harmonized_metadata_v2_NYU2.xlsx")
OUT = os.path.join(ROOT, "deliverables", "nyu2_impact.json")
os.makedirs(os.path.dirname(OUT), exist_ok=True)

ctx = build_context(Config())
atlas = set(ctx.tables["pseudobulks"]["sample_key"].astype(str))

xl = pd.ExcelFile(XLSX)
MM = xl.parse("Mutation_Matrix")
MM["Dataset"] = MM["Dataset"].astype(str).str.strip(); MM["Sample"] = MM["Sample"].astype(str).str.strip()
n2 = MM[MM["Dataset"] == "NYU-2"].copy()
n2["key"] = "NYU-2::" + n2["Sample"]
in_atlas = n2[n2["key"].isin(atlas)]
missing = sorted(set(n2["key"]) - atlas)
gcols = [c for c in MM.columns if c not in ("Dataset", "Sample", "key")]
X = in_atlas[gcols].apply(pd.to_numeric, errors="coerce").fillna(0)
new_pos = {g: int(X[g].sum()) for g in gcols if X[g].sum() > 0}
print("NYU-2 samples in file: %d | WITH atlas expression: %d | without: %d %s"
      % (len(n2), len(in_atlas), len(missing), missing))
print("new positives (atlas-present samples only): %d across %d genes" % (sum(new_pos.values()), len(new_pos)))

# ---- current per-mutation n_pos from the pipeline matrix ----
rows = list(csv.reader(open(os.path.join(HERE, "mutation_matrix_explicit_v2.tsv")), delimiter="\t"))
hdr = rows[0]
cols = [(i, c) for i, c in enumerate(hdr) if c.startswith(("mut_", "cyto_"))]
cur = {}
for i, c in cols:
    short = c.replace("mut_", "").replace("cyto_", "")
    n = sum(1 for r in rows[1:] if r[0] in atlas and r[i].strip() not in ("", "0", "0.0", "nan", "NA"))
    cur[short.upper()] = (c, n)

MINP = 8
tbl = []
for gene, add in sorted(new_pos.items(), key=lambda kv: -kv[1]):
    G = gene.upper().replace("-", "").replace("_", "")
    match = None
    for k in cur:
        if k.replace("-", "").replace("_", "") == G:
            match = k; break
    old = cur[match][1] if match else 0
    tbl.append({"gene": gene, "matched_field": cur[match][0] if match else None,
                "current_n_pos": old, "nyu2_adds": add, "new_n_pos": old + add,
                "was_trainable": old >= MINP, "now_trainable": (old + add) >= MINP})
newly = [t for t in tbl if t["now_trainable"] and not t["was_trainable"]]
print("\n%-12s %8s %8s %8s   %s" % ("gene", "current", "+NYU-2", "new", "status"))
for t in tbl:
    flag = "NEWLY TRAINABLE" if (t["now_trainable"] and not t["was_trainable"]) else ("trainable" if t["now_trainable"] else "still underpowered")
    print("%-12s %8d %8d %8d   %s" % (t["gene"], t["current_n_pos"], t["nyu2_adds"], t["new_n_pos"], flag))
print("\nmutations crossing the >=%d-positive threshold: %s" % (MINP, [t["gene"] for t in newly] or "none"))

json.dump({"n_file_samples": len(n2), "n_with_atlas_expression": len(in_atlas), "missing_from_atlas": missing,
           "total_new_positives": int(sum(new_pos.values())), "per_gene": tbl,
           "newly_trainable": [t["gene"] for t in newly],
           "note": "NYU-2 previously had ALL-ZERO mutation rows (silently wild-type). These positives were "
                   "counted as true negatives, which understates specificity and removes real positives "
                   "from training."}, open(OUT, "w"), indent=1)
print("\nwrote", OUT)
print("NYU2 IMPACT OK")
