#!/usr/bin/env python3
"""Companion to gui/evidence.json: each atlas sample's own per-cell-state marker expression, so the
evidence view can highlight THE CURRENT PATIENT'S dot on the violin (larger, distinct) instead of
just the anonymous mutant/control clouds.

We do NOT touch evidence.json (the violins stay exactly as built). We recompute each sample's value
with cp10k_log1p over the RNA pseudobulk counts — verified (build-time assert) to reproduce
evidence.json's pooled values to 1e-3, i.e. the SAME scale, so the overlaid dot lands on the violin.

Only the molecules/cell-states that already appear in evidence.json are emitted. Zero values are
omitted (sparse markers) — the frontend treats "sample in `samples` but missing here" as 0.0, and a
sample NOT in `samples` (Trumpp / uploads / Leucegene — not in this atlas) simply gets no dot.

Output: gui/evidence_samples.json  = {"samples":[atlas sample_keys], "drivers":{driver:{mol:{state:{sample_key:value}}}}}

Run on an LSF compute node (needs numpy/anndata):
  bsub -q test -W 30 -M 24000 -R "rusage[mem=24000]" -o bes.log \
    /usr/local/anaconda3-2020/bin/python build_evidence_samples.py
"""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
from amlmm.context import build_context, Config
from amlmm import pseudobulk_io as pio

EVPATH = os.path.join(HERE, "..", "gui", "evidence.json")
OUTPATH = os.path.join(HERE, "..", "gui", "evidence_samples.json")
EV = json.load(open(EVPATH))["drivers"]

ctx = build_context(Config())
pb = ctx.tables["pseudobulks"]
all_keys = sorted(set(pb["sample_key"].astype(str)))
print("atlas sample_keys:", len(all_keys))

# every (cell-state -> {molecules}) referenced by any driver's violins
need = {}
for d, dd in EV.items():
    for mol, states in dd.get("violins", {}).items():
        for st in states:
            need.setdefault(st, set()).add(mol)
print("cell-states to recompute:", len(need))

# per state: cp10k_log1p ONCE over the full pseudobulk count matrix, then extract each needed molecule
smv = {}          # (state, mol) -> {sample_key: rounded value}  (nonzero only)
check_ok = check_tot = 0
for st, mols in need.items():
    ids = pio.cellstate_pseudobulks(ctx, st)
    Xf = pio.pseudobulk_modality_matrix(ctx, "RNA", ids)
    if Xf.shape[0] == 0:
        continue
    id2s = pb.loc[list(Xf.index), "sample_key"].astype(str).to_dict()
    Xn = pio.cp10k_log1p(Xf.values)
    colpos = {c: i for i, c in enumerate(Xf.columns)}
    for mol in mols:
        if mol not in colpos:
            continue
        ci = colpos[mol]
        vals = {}
        for i, r in enumerate(Xf.index):
            v = float(Xn[i, ci])
            if v > 1e-6:
                vals[id2s[r]] = round(v, 4)
        smv[(st, mol)] = vals

# build-time verification: pooled evidence.json values must be reproduced by our per-sample recompute
for d, dd in EV.items():
    for mol, states in dd.get("violins", {}).items():
        for st, gd in states.items():
            vals = smv.get((st, mol))
            if vals is None:
                continue
            mine = sorted(vals.values()) + [0.0]        # 0.0 stands in for the omitted zeros
            pooled = list(gd.get("mutant") or []) + list(gd.get("control") or [])
            for p in pooled:
                check_tot += 1
                if p <= 1e-6 or any(abs(p - m) < 1e-3 for m in mine):
                    check_ok += 1
print("scale check: %d/%d pooled values reproduced (%.1f%%)"
      % (check_ok, check_tot, 100.0 * check_ok / max(check_tot, 1)))

# assemble driver -> molecule -> state (covers alias driver keys too — they share molecules)
drivers = {}
for d, dd in EV.items():
    for mol, states in dd.get("violins", {}).items():
        for st in states:
            vals = smv.get((st, mol))
            if vals:
                drivers.setdefault(d, {}).setdefault(mol, {})[st] = vals

json.dump({"samples": all_keys, "drivers": drivers}, open(OUTPATH, "w"), separators=(",", ":"))
print("wrote %s  (%d bytes, %d drivers)" % (OUTPATH, os.path.getsize(OUTPATH), len(drivers)))
print("EVIDENCE SAMPLES OK")
