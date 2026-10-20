#!/usr/bin/env python3
"""One row per held-out call, with both of its possible labels, for `rigor labels`.

  patient   the report's run folder
  cohort    sealed held-out scRNA, Trumpp/Waclawiczek, GSE281087
  model     multimodal (sealed set) or bulk caller (the other two)
  mutation  the category the call is for
  call      present | absent
  recorded  the report's true_label, blank when none is recorded
  inferred  pooled_heldout_eval.py's label from the cohort's driver list,
            blank when the cohort lists no drivers

Writes deliverables/heldout_calls.tsv. Then:

  rigor labels deliverables/heldout_calls.tsv --group cohort
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import pooled_heldout_variants as v

OUT = os.path.join(ROOT, "deliverables", "heldout_calls.tsv")
with open(OUT, "w", encoding="utf-8", newline="\n") as fh:
    fh.write("patient\tcohort\tmodel\tmutation\tcall\trecorded\tinferred\n")
    n = 0
    for p in v.load():
        for c in p["calls"]:
            inferred = ""
            if p["has_known"]:
                inferred = "present" if v.shipped_rule(c["name"], p["known"]) else "absent"
            fh.write("\t".join([p["id"], p["cohort"], p["model"], c["name"],
                                "present" if c["call"] else "absent", c["tl"] or "",
                                inferred]) + "\n")
            n += 1
print("wrote deliverables/heldout_calls.tsv (%d calls)" % n)
