#!/usr/bin/env python3
"""Build the standalone bake-off bundle: run bulk_external.harmonize() -> bundle_data.npz (variant-level labels)."""
import os, sys, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "pipeline"))
import numpy as np
import bulk_external as BE

OUT = r"C:\Users\krog5w\.gemini\antigravity\scratch\aml-bakeoff\bundle_data.npz"

t0 = time.time()
H = BE.harmonize()
genes = np.array(H["genes"], dtype=object)
drivers = np.array(BE.DRIVERS, dtype=object)
holdout = np.array(sorted(H["holdout"]), dtype=object)

pack = {"genes": genes, "drivers": drivers, "holdout": holdout}
for name, p in [("beataml", "ba"), ("leucegene", "lg"), ("sc", "sc")]:
    Xlin = H[name]["Xlin"]                                   # DataFrame samples x genes (linear)
    L = H[name]["L"].reindex(Xlin.index)[BE.DRIVERS]         # samples x drivers, aligned
    pack[p + "_samples"] = np.array([str(s) for s in Xlin.index], dtype=object)
    pack[p + "_X"] = Xlin.values.astype(np.float32)
    pack[p + "_L"] = L.values.astype(np.float32)
    print("  %-10s X=%s L=%s  labeled/pos per-driver saved" % (name, Xlin.shape, L.shape))

np.savez_compressed(OUT, **pack)
mb = os.path.getsize(OUT) / 1e6
print("wrote %s  (%.1f MB, %d genes, %d categories, %.0fs)" % (OUT, mb, len(genes), len(drivers), time.time() - t0))
# quick sanity: how many categories clear >=6 positives in BeatAML
baL = H["beataml"]["L"]
elig = [d for d in BE.DRIVERS if int((baL[d] == 1).sum()) >= 6 and int((baL[d] == 0).sum()) >= 6]
print("BeatAML CV-eligible (>=6 pos & >=6 neg): %d categories" % len(elig))
