"""Model B — the same response model, applied one malignant cell state at a time.

BeatAML's own analysis concluded that differentiation state broadly determines ex-vivo drug response.
A bulk average therefore hides the thing that matters most: a 5% primitive subpopulation predicted
resistant is more clinically relevant than the 80% monocytic bulk predicted sensitive, and the bulk
number cannot tell you the difference.

Because Model A consumes one expression vector in the shared gene space, it can be handed a single
cell state's pseudobulk exactly as it would be handed a whole sample. That gives, per drug:

    sens_k            predicted response of cell state k
    sens_bulk         prediction from the whole-sample bulk-equivalent (what a bulk assay would see)
    sens_weighted     sum_k pi_k * sens_k -- the abundance-weighted mean
    coverage_*        fraction of the (presumed) malignant / LSC-like compartment predicted sensitive
    worst_state       the most resistant state above an abundance floor: the escape candidate
    dispersion        spread of response across states
    bulk_vs_sc        sens_bulk - sens_weighted; large values mean the bulk view is misleading here

**What "malignant" means here.** We do not have per-cell somatic genotypes, so we do not claim to
identify malignant cells. We use the *blast/progenitor compartment* -- the states in which AML blasts
reside -- as the candidate compartment, and label it as such everywhere it is reported. Lymphoid and
stromal states are treated as presumed-normal bystanders and excluded from coverage.
"""
from __future__ import annotations
import re
import numpy as np
import pandas as pd

# 89 atlas states -> lineage group. Kept identical to build_state_signatures.GROUPS so the signature
# scores and the compartment definitions cannot drift apart.
GROUPS = [
    ("stroma",       r"^(MSC|Osteoblast|Stromal)"),
    ("T_NK",         r"^(CD4 |CD8 |MAIT|NK-|T CD)"),
    ("B_lineage",    r"^(Pro-B|pre-B|Transitional-B|B Memory|Plasma Cell|BMCP)"),
    ("DC",           r"^(cDC|pDC|ASDC|pre-DC)"),
    ("erythroid",    r"^(ERP-|Erythroblast|MEP-Eryth)"),
    ("MEP_Mk",       r"^(MEP-|MKP-|MK-Platelet|MPP-MEP)"),
    ("monocytic",    r"^(Classical-Mono|Intermediate Mono|Mono-|Non-Classical Mono|Mac$|MDP-)"),
    ("granulocytic", r"^(cMOP|preNeu|immNeu|Myeloid intermediate)"),
    ("GMP",          r"^(MultiLin-GMP)"),
    ("LMPP_CLP",     r"^(LMPP|CLP)"),
    ("HSC_MPP",      r"^(HSC-|MPP-|Multilin-)"),
]
# the compartment AML blasts occupy; coverage is reported over this, never over all cells
BLAST_GROUPS = {"HSC_MPP", "LMPP_CLP", "GMP", "granulocytic", "monocytic", "MEP_Mk", "erythroid", "DC"}
# the primitive compartment used as the LSC-like proxy
LSC_GROUPS = {"HSC_MPP", "LMPP_CLP", "GMP"}
LSC_STATES = re.compile(r"^(HSC-|MPP-|LMPP|MultiLin-GMP|Multilin-)")

MIN_CELLS = 20          # a state with fewer cells has an unusable pseudobulk
MIN_FRAC = 0.01         # states below 1% of the sample are excluded from worst-state / escape calls


def group_of(state):
    for name, pat in GROUPS:
        if re.match(pat, str(state)):
            return name
    return None


def cp10k(counts):
    """(rows x gene) raw counts -> LINEAR counts-per-10k. Linear because that is the scale the shared
    bundle stores (`ba_X`/`sc_X`); the model applies its own log2 and z-scoring downstream."""
    v = np.asarray(counts, dtype=np.float64)
    tot = v.sum(1, keepdims=True); tot[tot == 0] = 1.0
    return v / tot * 1e4


def state_pseudobulks(X_counts, genes, labels, min_cells=MIN_CELLS):
    """Sum RAW counts per assigned cell state. Raw first, normalise after -- summing pre-normalised
    values would weight a 30-cell state the same as a 3,000-cell one."""
    lab = pd.Series(np.asarray(labels))
    out, n = {}, {}
    X = X_counts
    for st, idx in lab.groupby(lab).groups.items():
        ii = np.asarray(idx)
        if len(ii) < min_cells:
            continue
        sub = X[ii]
        v = np.asarray(sub.sum(0)).ravel() if hasattr(sub, "sum") else np.asarray(sub).sum(0)
        out[st] = v
        n[st] = int(len(ii))
    if not out:
        return pd.DataFrame(columns=list(genes)), pd.Series(dtype=int)
    return pd.DataFrame(out, index=list(genes)).T, pd.Series(n)


class StateResponse:
    """Wraps a trained Model A and scores each cell state of one patient."""

    def __init__(self, model, sym2ens=None):
        self.mod = model
        self.sym2ens = sym2ens or getattr(model, "sym2ens", {})

    def _z(self, expr_lin_df, ref="sc"):
        """Align an arbitrary (rows x genes) frame to the model's gene space and z-score it against a
        MATCHED cohort reference.

        Matched matters: single-cell bulk-equivalents and BeatAML bulk sit on different platforms with
        different gene-wise means, and z-scoring a single-cell sample against BeatAML's means pushes it
        off the end of the training distribution for every gene at once. Falls back to whatever
        reference exists if the requested one was never built."""
        fs = self.mod.fs
        # Column->column index map, cached per input gene list. Building this as a DataFrame and
        # assigning 14,237 columns is ~100x slower than one numpy gather, and this runs once per cell
        # state per patient.
        cols_ = expr_lin_df.columns
        key = (len(cols_), str(cols_[0]), str(cols_[-1]))     # id() would be reused after GC
        cache = getattr(self, "_align_cache", None)
        if cache is None:
            cache = self._align_cache = {}
        plan = cache.get(key)
        if plan is None:
            gset = set(fs.genes)
            gpos = {g: i for i, g in enumerate(fs.genes)}
            src, dst = [], []
            seen = set()
            for j, g in enumerate(expr_lin_df.columns):
                g = str(g)
                e = g if g in gset else self.sym2ens.get(g)
                if e in gpos and e not in seen:
                    seen.add(e); src.append(j); dst.append(gpos[e])
            plan = cache[key] = (np.asarray(src), np.asarray(dst))
        src, dst = plan
        A = np.zeros((len(expr_lin_df), len(fs.genes)), dtype=np.float64)
        if len(src):
            A[:, dst] = np.asarray(expr_lin_df.values, dtype=np.float64)[:, src]
        use = ref if ref in fs.ref else ("beataml" if "beataml" in fs.ref else next(iter(fs.ref)))
        return fs.z(A, use), len(src)

    def predict(self, state_counts, n_cells, mut=None, meta=None, drugs=None,
                min_frac=MIN_FRAC, prob_cut=0.5, state_ref="sc_state", bulk_ref="sc_sample"):
        """state_counts: (state x gene) RAW summed counts. n_cells: cells per state.

        Returns {'per_state': {state: {drug: sens}}, 'per_drug': {drug: {...aggregates...}}}.
        A cell state is percentiled against OTHER CELL STATES (`sc_state`) and the whole sample against
        OTHER WHOLE SAMPLES (`sc_sample`) -- mixing the two would make every state look like an outlier.
        """
        if not len(state_counts):
            return {"per_state": {}, "per_drug": {}, "note": "no cell state passed the cell-count floor"}
        expr = pd.DataFrame(cp10k(state_counts.values), index=state_counts.index,
                            columns=state_counts.columns)
        Z, n_shared = self._z(expr)
        states = list(state_counts.index)
        grp = {s: group_of(s) for s in states}
        n_cells = pd.Series(n_cells).reindex(states).fillna(0)
        blast = np.array([grp[s] in BLAST_GROUPS for s in states])
        lsc = np.array([bool(LSC_STATES.match(str(s))) for s in states])
        frac_all = (n_cells / max(1.0, n_cells.sum())).values
        bl_tot = n_cells.values[blast].sum()
        lsc_tot = n_cells.values[lsc].sum()

        # whole-sample bulk-equivalent: sum raw counts over ALL states, THEN normalise
        bulk_lin = cp10k(state_counts.values.sum(0, keepdims=True))
        Zb, _ = self._z(pd.DataFrame(bulk_lin, index=["bulk"], columns=state_counts.columns))

        # the patient feature matrices do not depend on the drug -- build them once
        Xp, _, _ = self.mod.fs.transform(Z, meta=None, mut=None, blocks=self.mod.blocks)
        Xb, _, _ = self.mod.fs.transform(Zb, meta=meta, mut=mut, blocks=self.mod.blocks)
        allrows = np.arange(len(states))

        per_state, per_drug = {}, {}
        drugs = drugs or list(self.mod.drug_models)
        for drug in drugs:
            dm = self.mod.drug_models.get(drug)
            gm = self.mod.group_models.get(self.mod.drug_group.get(drug))
            if dm is None or gm is None:
                continue
            T, _ = self.mod.fs.target_features(Z, _targets(drug), self.sym2ens)
            s = self.mod._shared_pred(gm, Xp, allrows, T) + dm["w"] * dm["est"].predict(Xp)

            Tb, _ = self.mod.fs.target_features(Zb, _targets(drug), self.sym2ens)
            sb = float(self.mod._shared_pred(gm, Xb, np.array([0]), Tb)[0] + dm["w"] * dm["est"].predict(Xb)[0])

            p = np.array([self.mod.calibrated(drug, x, state_ref) or np.nan for x in s])
            sens_state = {st: round(float(x), 4) for st, x in zip(states, s)}
            w_all = float(np.nansum(frac_all * s))
            hit = np.nan_to_num(p, nan=0.0) >= prob_cut
            cov_b = (float(n_cells.values[blast & hit].sum() / bl_tot) if bl_tot else None)
            cov_l = (float(n_cells.values[lsc & hit].sum() / lsc_tot) if lsc_tot else None)
            elig = blast & (frac_all >= min_frac)
            worst = (states[int(np.nanargmin(np.where(elig, s, np.nan)))] if elig.any() else None)
            per_drug[drug] = {
                "sens_bulk": round(sb, 4),
                "prob_sensitive_bulk": self.mod.calibrated(drug, sb, bulk_ref),
                "percentile_bulk": self.mod.score_percentile(drug, sb, bulk_ref),
                "state_probabilities": {st: (None if not np.isfinite(x) else round(float(x), 3))
                                        for st, x in zip(states, p)},
                "sens_weighted": round(w_all, 4),
                "bulk_vs_sc": round(sb - w_all, 4),
                "coverage_blast": None if cov_b is None else round(cov_b, 3),
                "coverage_LSC_like": None if cov_l is None else round(cov_l, 3),
                "worst_state": worst,
                "worst_state_sens": None if worst is None else round(float(s[states.index(worst)]), 4),
                "dispersion": round(float(np.nanstd(s[blast])), 4) if blast.any() else None,
                "n_states_scored": int(np.isfinite(s).sum()),
            }
            per_state[drug] = sens_state
        return {"per_state": per_state, "per_drug": per_drug,
                "states": [{"state": s, "group": grp[s], "n_cells": int(n_cells[s]),
                            "fraction": round(float(f), 4), "blast_compartment": bool(b),
                            "LSC_like": bool(l)}
                           for s, f, b, l in zip(states, frac_all, blast, lsc)],
                "genes_shared": int(n_shared),
                "compartment_note": ("blast compartment = the states AML blasts occupy; no per-cell "
                                     "somatic genotype is used, so this is a compartment, not a "
                                     "malignant-cell call")}


def _targets(drug):
    from . import targets as TG
    return TG.get(drug)["targets"]
