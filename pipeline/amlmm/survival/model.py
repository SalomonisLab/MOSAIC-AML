"""The survival layer: per-block Cox models, stacked, with the clinical baseline kept separate on purpose.

Structure mirrors COMPASS-AML, for the same reason: the feature blocks differ by orders of magnitude in
width (100 RNA PCs versus 4 age/ELN columns) and one shared penalty over the concatenation lets the wide
block dominate whether or not it deserves to.

    risk_i = stack( cox_rna(x_i), cox_state(x_i), cox_mut(x_i), cox_clin(x_i) )

Each block gets its own ridge-penalised Cox; their linear predictors become the four covariates of a
small stacking Cox fitted on **inner out-of-fold** predictions, so the stacker never sees a block's
in-sample optimism. If the stack does not beat the best single block on those same held-out rows, the
weights collapse onto that block — fusion can never do worse than its best part.

The question this layer exists to answer is not "what is the C-index" but **"does anything molecular add
to age and ELN 2017"**, so `MOLECULAR_BLOCKS` is kept distinct from the clinical block and the headline
comparison is always incremental.
"""
from __future__ import annotations
import numpy as np

from .coxph import CoxPH, c_index

MOLECULAR_BLOCKS = ("rna", "state", "mut")
ALL_BLOCKS = ("rna", "state", "mut", "clin")


class SurvivalModel:
    VERSION = "1.0"

    def __init__(self, alphas=None, inner_folds=3):
        self.alphas = alphas or {"rna": 0.30, "state": 0.10, "mut": 0.20, "clin": 0.05}
        self.inner_folds = int(inner_folds)
        self.blocks = []
        self.models = {}          # block -> CoxPH
        self.stack = None         # CoxPH over the block linear predictors
        self.stack_blocks = []
        self.floored_to = None
        self.dead_blocks = []
        self.risk_ref = None      # sorted training risks, for percentiles
        self.meta = {}

    # ------------------------------------------------------------ fit ----
    def fit(self, blocks, time, event, groups):
        """blocks: {name: X}. Inner grouped CV -> OOF block risks -> stacking Cox, with the floor guard."""
        from sklearn.model_selection import GroupKFold
        self.blocks = [b for b in ALL_BLOCKS if b in blocks]
        t, e, g = np.asarray(time, float), np.asarray(event, int), np.asarray(groups)
        n = len(t)

        oof = {b: np.full(n, np.nan) for b in self.blocks}
        nsp = min(self.inner_folds, len(np.unique(g)))
        if nsp >= 2:
            for tri, tei in GroupKFold(n_splits=nsp).split(np.zeros(n), groups=g):
                for b in self.blocks:
                    m = CoxPH(alpha=self.alphas.get(b, 0.1)).fit(blocks[b][tri], t[tri], e[tri])
                    oof[b][tei] = m.risk(blocks[b][tei])
        # A block that failed to fit leaves an all-NaN OOF column; keeping it in the "all blocks
        # scored" mask empties the mask and silently collapses the fusion onto an arbitrary block.
        # Drop such blocks outright and record it.
        dead = [b for b in self.blocks if not np.isfinite(oof[b]).any()]
        if dead:
            self.blocks = [b for b in self.blocks if b not in dead]
        self.dead_blocks = dead
        if not self.blocks:
            raise RuntimeError("every block failed to fit")
        ok = np.all([np.isfinite(oof[b]) for b in self.blocks], axis=0)

        singles = {b: c_index(t[ok], e[ok], oof[b][ok]) for b in self.blocks}
        best_b = max(singles, key=lambda b: singles[b])
        A = np.column_stack([oof[b][ok] for b in self.blocks])
        stacked_c = float("nan")
        if ok.sum() >= 40:
            st = CoxPH(alpha=0.05).fit(A, t[ok], e[ok])
            stacked_c = c_index(t[ok], e[ok], st.risk(A))

        if not (stacked_c == stacked_c) or stacked_c <= singles[best_b]:
            self.stack_blocks = [best_b]
            self.floored_to = best_b
            self.stack = None
        else:
            self.stack_blocks = list(self.blocks)
            self.floored_to = None
            self.stack = CoxPH(alpha=0.05).fit(A, t[ok], e[ok])

        for b in self.blocks:                                   # deployed per-block models: all data
            self.models[b] = CoxPH(alpha=self.alphas.get(b, 0.1)).fit(blocks[b], t, e)
        self.meta = {"blocks": self.blocks, "inner_oof_cindex": {k: round(v, 4) for k, v in singles.items()},
                     "stacked_inner_cindex": None if stacked_c != stacked_c else round(stacked_c, 4),
                     "floored_to": self.floored_to, "dead_blocks": self.dead_blocks,
                     "n": int(n), "events": int(e.sum())}
        self._set_reference(blocks, t, e)
        return self

    def _set_reference(self, blocks, t, e):
        r = self.risk(blocks)
        self.risk_ref = np.sort(r)
        # a baseline hazard on the DEPLOYED risk scale, so survival curves match the reported risk
        self._final = CoxPH(alpha=1e-6).fit(r.reshape(-1, 1), t, e)

    # -------------------------------------------------------- predict ----
    def risk(self, blocks):
        parts = [self.models[b].risk(blocks[b]) for b in self.stack_blocks]
        if self.stack is None:
            return parts[0]
        return self.stack.risk(np.column_stack(parts))

    def percentile(self, r):
        if self.risk_ref is None:
            return None
        r = np.atleast_1d(r)
        return np.searchsorted(self.risk_ref, r, side="right") / len(self.risk_ref)

    def survival(self, blocks, times):
        return self._final.survival(self.risk(blocks).reshape(-1, 1), times)

    def median_survival(self, blocks):
        return self._final.median_survival(self.risk(blocks).reshape(-1, 1))
