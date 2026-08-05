"""Model C — is the drug's target pathway actually active in this patient?

Deliberately kept OUT of the empirical response model. Model A learns "patients who look like this
responded ex vivo"; Model C asks "is there a mechanistic reason this compound should work here?".
Keeping them separate is the whole point: agreement between an empirical and a mechanistic line of
evidence is much stronger than either alone, and disagreement is informative rather than something to
average away. Fusing them into one score would destroy exactly the signal that makes the pair useful.

Evidence assembled per (patient, drug):

  target        expression percentile of the compound's own annotated targets
  pathway       a short curated transcriptional readout of the pathway the compound is aimed at --
                output genes, not the pathway members, since a kinase's mRNA says little about its
                activity while its transcriptional output says a great deal
  dependency    the anti-apoptotic balance BCL2 / (BCL2 + MCL1 + BCL2L1 + BCL2A1), the best-evidenced
                mechanistic determinant of BH3-mimetic response
  genetic       an activating lesion in the target pathway, taken from the OBSERVED genotype where
                available and otherwise from the calibrated mutation caller (flagged as such)
  resistance    measurable proxies for this compound's curated bypass routes

Percentiles are against the BeatAML training cohort, so "high" means high among AML, not high in the
abstract.
"""
from __future__ import annotations
import numpy as np
import pandas as pd

from . import targets as TG

# short transcriptional OUTPUT readouts, per target-pathway family group
PATHWAY_READOUT = {
    "MAPK":              ["DUSP4", "DUSP6", "SPRY2", "SPRY4", "ETV4", "ETV5", "PHLDA1"],
    "FLT3":              ["DUSP6", "SPRY2", "PIM1", "MYC", "CISH", "SOCS2"],
    "RTK":               ["DUSP6", "SPRY2", "SPRY4", "ETV5", "PHLDA1"],
    "JAK_STAT":          ["SOCS1", "SOCS2", "SOCS3", "CISH", "PIM1", "BCL2L1"],
    "PI3K_AKT_mTOR":     ["RPS6KB1", "EIF4EBP1", "SLC7A5", "SLC2A1", "HK2", "TXNIP"],
    "NF-kB":             ["NFKBIA", "TNFAIP3", "BIRC3", "CCL2", "CXCL8", "RELB"],
    "cell_cycle":        ["MKI67", "TOP2A", "CCNB1", "CDK1", "BUB1", "PLK1", "AURKB", "E2F1"],
    "apoptosis":         ["BCL2", "MCL1", "BCL2L1", "BCL2A1", "BAX", "BAK1", "BCL2L11", "BIK"],
    "epigenetic":        ["MYC", "MYB", "HEXIM1", "CDK6", "BCL2"],
    "proteostasis":      ["XBP1", "HSPA5", "DDIT3", "HSPA1A", "HSPH1"],
    "immune_signalling": ["BTK", "SYK", "LYN", "CD79A", "CD79B", "BLNK"],
    "developmental":     ["HES1", "HEY1", "AXIN2", "LEF1", "SMAD7", "SERPINE1"],
    "metabolic":         ["IDH1", "IDH2", "SLC7A11", "GPX4", "FDX1"],
    "chemotherapy":      ["DCK", "SLC29A1", "CDA", "SAMHD1", "RRM1", "RRM2"],
    "stress_MAPK":       ["JUN", "FOS", "ATF3", "DUSP1", "GADD45B"],
}

# curated resistance string -> a measurable proxy. `direction` is the direction that means RESISTANT.
RESISTANCE_PROXY = {
    "MCL1 upregulation":            ("expr", ["MCL1"], "high"),
    "BCL2L1 (BCL-xL)":              ("expr", ["BCL2L1"], "high"),
    "monocytic differentiation":    ("state", ["monocytic"], "high"),
    "TP53 loss":                    ("mut", ["TP53"], "present"),
    "TP53 mutation/loss (absolute)": ("mut", ["TP53"], "present"),
    "BAX loss":                     ("expr", ["BAX"], "low"),
    "RB1 loss":                     ("expr", ["RB1"], "low"),
    "CCNE1 amplification":          ("expr", ["CCNE1"], "high"),
    "DCK loss":                     ("expr", ["DCK"], "low"),
    "SLC29A1 (hENT1) loss":         ("expr", ["SLC29A1"], "low"),
    "CDA overexpression":           ("expr", ["CDA"], "high"),
    "SAMHD1 high":                  ("expr", ["SAMHD1"], "high"),
    "CRBN loss/mutation":           ("expr", ["CRBN"], "low"),
    "RAS/MAPK activation":          ("mut", ["NRAS", "KRAS", "PTPN11", "NF1"], "present"),
    "AXL/GAS6 bypass":              ("expr", ["AXL", "GAS6"], "high"),
    "MET amplification":            ("expr", ["MET"], "high"),
    "ERBB3 upregulation":           ("expr", ["ERBB3"], "high"),
    "MCL1 re-expression":           ("expr", ["MCL1"], "high"),
    "eIF4E amplification":          ("expr", ["EIF4E"], "high"),
    "SGK1 compensation":            ("expr", ["SGK1"], "high"),
    "XIAP":                         ("expr", ["XIAP"], "high"),
    "IAP redundancy (XIAP)":        ("expr", ["XIAP"], "high"),
    "WNT activation":               ("expr", ["AXIN2", "LEF1"], "high"),
    "proteasome subunit upregulation": ("expr", ["PSMB5", "PSMA1", "PSMB1"], "high"),
    "HSPA1A/HSF1 heat-shock response": ("expr", ["HSPA1A", "HSF1"], "high"),
    "glycolytic shift":             ("expr", ["HK2", "SLC2A1", "LDHA"], "high"),
    "low FDX1":                     ("expr", ["FDX1"], "low"),
    "antioxidant capacity (GPX4/NRF2)": ("expr", ["GPX4", "NFE2L2"], "high"),
}

# variant-level driver categories in the bundle that activate each pathway
GENETIC_ACTIVATION = {
    "FLT3": ["FLT3_ITD", "FLT3_TKD", "FLT3_N676", "FLT3_other"],
    "MAPK": ["NRAS", "KRAS", "PTPN11", "NF1", "BRAF", "CBL"],
    "RTK":  ["KIT", "FLT3_ITD", "NRAS", "KRAS", "CBL"],
    "JAK_STAT": ["JAK2", "MPL", "CALR", "SH2B3"],
    "metabolic": ["IDH1", "IDH2"],
    "epigenetic": ["DNMT3A", "TET2", "ASXL1", "EZH2", "KMT2A"],
    "apoptosis": ["TP53"],
    "immune_signalling": [],
    "cell_cycle": ["TP53"],
}


class MechanismModel:
    """Percentile-references are built once from the BeatAML training matrix."""

    def __init__(self, genes, X_ref_lin, sym2ens=None, state_scores_ref=None):
        self.genes = list(genes)
        self.gidx = {g: i for i, g in enumerate(self.genes)}
        self.sym2ens = sym2ens or {}
        L = np.log2(np.clip(np.asarray(X_ref_lin, float), 0, None) + 1.0)
        self.ref_sorted = np.sort(L, axis=0)              # per-gene training distribution
        self.state_ref = state_scores_ref                 # DataFrame (samples x state score) or None

    # -------------------------------------------------------------- utils ----
    def _cols(self, symbols):
        out = []
        for s in symbols:
            e = s if s in self.gidx else self.sym2ens.get(s)
            if e in self.gidx:
                out.append(self.gidx[e])
        return out

    def _pct(self, x_lin, cols):
        """Percentile of this sample's expression against the training cohort, per gene, averaged."""
        if not cols:
            return None
        L = np.log2(np.clip(np.asarray(x_lin, float), 0, None) + 1.0)
        ps = []
        n = self.ref_sorted.shape[0]
        for c in cols:
            ps.append(np.searchsorted(self.ref_sorted[:, c], L[c], side="right") / n)
        return float(np.mean(ps))

    # ------------------------------------------------------------ scoring ----
    def evaluate(self, drug, x_lin, mutations=None, state_scores=None):
        """mutations: {category -> 0/1/probability} plus optional '__observed__' set naming which
        categories came from a real genotype rather than the caller."""
        ann = TG.get(drug)
        fam = ann["family_group"]
        ev = {"inhibitor": drug, "family_group": fam, "targets": ann["targets"]}

        ev["target_expression_pct"] = self._pct(x_lin, self._cols(ann["targets"]))
        ev["pathway_readout_pct"] = self._pct(x_lin, self._cols(PATHWAY_READOUT.get(fam, [])))

        # anti-apoptotic balance -- reported for every drug, decisive only for BH3 mimetics
        bcl2 = self._cols(["BCL2"]); comp = self._cols(["MCL1", "BCL2L1", "BCL2A1"])
        if bcl2 and comp:
            L = np.log2(np.clip(np.asarray(x_lin, float), 0, None) + 1.0)
            num = float(L[bcl2[0]]); den = num + float(np.sum(L[comp]))
            ev["bcl2_dependency"] = round(num / den, 4) if den > 0 else None
        else:
            ev["bcl2_dependency"] = None

        # genetic activation of the target pathway
        gen, obs = [], set((mutations or {}).get("__observed__", []) or [])
        for cat in GENETIC_ACTIVATION.get(fam, []):
            for k, v in (mutations or {}).items():
                if k == "__observed__" or not isinstance(v, (int, float)):
                    continue
                if k.upper().startswith(cat.upper()) and v >= 0.5:
                    gen.append({"lesion": k, "value": float(v),
                                "source": "observed" if k in obs else "predicted"})
        ev["genetic_activation"] = gen

        # resistance proxies
        res = []
        for r in ann["resistance"]:
            spec = RESISTANCE_PROXY.get(r)
            if spec is None:
                res.append({"mechanism": r, "measurable": False})
                continue
            kind, keys, direction = spec
            if kind == "expr":
                p = self._pct(x_lin, self._cols(keys))
                if p is None:
                    res.append({"mechanism": r, "measurable": False}); continue
                hit = (p >= 0.75) if direction == "high" else (p <= 0.25)
                res.append({"mechanism": r, "measurable": True, "proxy": "+".join(keys),
                            "percentile": round(p, 3), "direction": direction, "flagged": bool(hit)})
            elif kind == "state":
                v = None if state_scores is None else state_scores.get(keys[0])
                p = None if (v is None or self.state_ref is None or keys[0] not in self.state_ref)\
                    else float((self.state_ref[keys[0]].values <= v).mean())
                res.append({"mechanism": r, "measurable": p is not None, "proxy": "state:" + keys[0],
                            "percentile": None if p is None else round(p, 3),
                            "flagged": bool(p is not None and p >= 0.75)})
            else:                                                     # mutation-based
                hit = any(float(v) >= 0.5 for k, v in (mutations or {}).items()
                          if k != "__observed__" and isinstance(v, (int, float))
                          and any(k.upper().startswith(g.upper()) for g in keys))
                res.append({"mechanism": r, "measurable": mutations is not None,
                            "proxy": "mut:" + "/".join(keys), "flagged": bool(hit)})
        ev["resistance"] = res
        ev["resistance_flags"] = int(sum(1 for r in res if r.get("flagged")))

        ev["mechanistic_score"] = self._combine(ev, fam)
        return ev

    @staticmethod
    def _combine(ev, fam):
        """A transparent, deliberately simple 0-1 roll-up. Anything more elaborate would imply a
        precision the underlying evidence does not have, and the itemised evidence is what the report
        actually shows the reader."""
        parts, wts = [], []
        if ev.get("target_expression_pct") is not None:
            parts.append(ev["target_expression_pct"]); wts.append(1.0)
        if ev.get("pathway_readout_pct") is not None:
            parts.append(ev["pathway_readout_pct"]); wts.append(1.0)
        if fam == "apoptosis" and ev.get("bcl2_dependency") is not None:
            parts.append(min(1.0, max(0.0, (ev["bcl2_dependency"] - 0.15) / 0.25))); wts.append(2.0)
        if ev.get("genetic_activation"):
            parts.append(1.0); wts.append(1.5)
        if not parts:
            return None
        base = float(np.average(parts, weights=wts))
        return round(max(0.0, base - 0.10 * ev.get("resistance_flags", 0)), 4)
