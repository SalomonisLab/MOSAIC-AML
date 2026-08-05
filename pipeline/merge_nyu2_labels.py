#!/usr/bin/env python3
"""Merge the new NYU-2 mutation calls into mutation_matrix_explicit_v2.tsv.

NYU-2's 28 samples previously had ALL-ZERO mutation rows — every gene read as wild-type, so 136 real
positives (in the 25 samples that have atlas expression) were being used as TRUE NEGATIVES. That is
label noise in training AND understated specificity in validation.

Only NYU-2 rows are touched; every other cohort's labels are byte-identical. Writes
mutation_matrix_explicit_v3.tsv (v2 kept intact) + a diff report.

  /usr/local/anaconda3-2020/bin/python merge_nyu2_labels.py [xlsx]
"""
import os, sys, csv, json
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
XLSX = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "labels", "AML_harmonized_metadata_v2_NYU2.xlsx")
SRC = os.path.join(HERE, "mutation_matrix_explicit_v2.tsv")
DST = os.path.join(HERE, "mutation_matrix_explicit_v3.tsv")
import pandas as pd

MM = pd.ExcelFile(XLSX).parse("Mutation_Matrix")
MM["Dataset"] = MM["Dataset"].astype(str).str.strip(); MM["Sample"] = MM["Sample"].astype(str).str.strip()
n2 = MM[MM["Dataset"] == "NYU-2"].copy()
n2["key"] = "NYU-2::" + n2["Sample"]
gcols = [c for c in MM.columns if c not in ("Dataset", "Sample", "key")]

rows = list(csv.reader(open(SRC), delimiter="\t"))
hdr = rows[0]
# map a matrix column (mut_X / cyto_X) to the xlsx gene column, normalizing punctuation
def norm(s): return str(s).upper().replace("-", "").replace("_", "").replace("(", "").replace(")", "").replace(";", "")
xls_by_norm = {norm(g): g for g in gcols}
col_map = {}
for i, c in enumerate(hdr):
    if not c.startswith(("mut_", "cyto_")):
        continue
    short = c.replace("mut_", "").replace("cyto_", "")
    g = xls_by_norm.get(norm(short))
    if g is not None:
        col_map[i] = g
print("matrix columns matched to the new file: %d" % len(col_map))
unmatched = [c for c in hdr if c.startswith(("mut_", "cyto_")) and norm(c.replace("mut_", "").replace("cyto_", "")) not in xls_by_norm]
print("  unmatched matrix columns (left unchanged): %s" % (unmatched or "none"))

new_by_key = {r["key"]: r for _, r in n2.iterrows()}
changed, flips = 0, {}
out = [hdr]
for r in rows[1:]:
    key = r[0]
    if key in new_by_key:
        src = new_by_key[key]; touched = False
        for i, g in col_map.items():
            v = src[g]
            try:
                v = int(float(v))
            except Exception:
                continue
            if v == 1 and r[i].strip() in ("", "0", "0.0", "nan", "NA"):
                r[i] = "1"; touched = True
                flips[hdr[i]] = flips.get(hdr[i], 0) + 1
        if touched:
            changed += 1
    out.append(r)

with open(DST, "w", newline="") as fh:
    csv.writer(fh, delimiter="\t").writerows(out)
print("\nrows changed: %d | total 0->1 flips: %d" % (changed, sum(flips.values())))
for k, v in sorted(flips.items(), key=lambda kv: -kv[1]):
    print("   %-22s +%d" % (k, v))
json.dump({"source": os.path.basename(XLSX), "rows_changed": changed, "flips": flips,
           "total_flips": sum(flips.values()), "output": os.path.basename(DST)},
          open(os.path.join(ROOT, "deliverables", "nyu2_merge_report.json"), "w"), indent=1)
print("\nwrote", DST)
print("MERGE NYU2 OK")
