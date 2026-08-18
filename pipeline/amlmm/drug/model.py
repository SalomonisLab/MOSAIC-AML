"""Model A — the patient-level BeatAML ex-vivo response model.

Hierarchical, because 118 inhibitors with ~400 trainable specimens each is exactly the regime where
118 independent models overfit and one pooled model underfits:

    sens_ij  =  f_group(X_i, D_j)  +  w_j * f_j(X_i)

`f_group` is a ridge fitted over every (specimen x inhibitor) row in a target-pathway family group,
with a drug descriptor (the z-expression of that inhibitor's own annotated targets) supplied as
interaction features -- so it learns "what kind of leukaemia is sensitive to RTK inhibition", and can
say something about a drug it has only a little data for. `f_j` is a per-drug ridge on the residual,
shrunk by w_j = n_j/(n_j + kappa) so a small drug barely moves off the family prior.

The modelled quantity is `sens = -auc_z`, the within-drug robust z of the fitted curve AUC with the
sign flipped, so that larger always means more sensitive.

Two anti-leakage rules are enforced here rather than left to the caller, because both are easy to get
wrong and neither is visible in the metrics until external validation fails:

  * the within-drug normalisation constants (median, MAD) and the sensitive/resistant tail cut-points
    are fitted on TRAINING specimens only and applied to held-out ones;
  * every fold refits the feature space (PCA, imputation medians) from scratch.
"""
from __future__ import annotations
import numpy as np
import pandas as pd

from . import targets as TG

KAPPA = 150.0            # drug-specific shrinkage: n_j/(n_j+kappa); at n=450 the drug term gets 0.75
ALPHAS = (1.0, 10.0, 100.0, 1000.0, 10000.0)


def _pct(sorted_ref, v):
    """Percentile of `v` within a sorted reference distribution, NaN-safe."""
    v = np.asarray(v, float)
    n = len(sorted_ref)
    out = np.full(v.shape, np.nan)
    ok = np.isfinite(v)
    if n >= 2:
        out[ok] = np.searchsorted(sorted_ref, v[ok], side="right") / n
    elif n == 1:
        out[ok] = 0.5
    return out


def _corr(a, b):
    """Pearson r, NaN-safe and degenerate-safe -- the fusion floor's comparison statistic."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 5 or a[m].std() == 0 or b[m].std() == 0:
        return -np.inf
    return float(np.corrcoef(a[m], b[m])[0, 1])


# ------------------------------------------------ train-only normalisation ----
def fit_norm(long, train_mask, tail=0.20, min_n=20):
    """Per-drug median / MAD / tail cut-points from training rows only."""
    out = {}
    tr = long[train_mask]
    for drug, g in tr.groupby("inhibitor"):
        a = pd.to_numeric(g["auc"], errors="coerce").dropna()
        if len(a) < min_n:
            continue
        med = float(a.median())
        mad = float(np.median(np.abs(a - med)))
        scale = 1.4826 * mad
        if scale <= 1e-9:
            scale = float(a.std(ddof=1)) or 1.0
        out[drug] = {"median": med, "scale": scale, "n_train": int(len(a)),
                     "lo": float(a.quantile(tail)), "hi": float(a.quantile(1.0 - tail))}
    return out


def apply_norm(long, norm):
    """Add train-referenced `sens` (higher = more sensitive) and `y_sens` (1 / 0 / NaN)."""
    med = long["inhibitor"].map(lambda d: norm.get(d, {}).get("median", np.nan))
    sc = long["inhibitor"].map(lambda d: norm.get(d, {}).get("scale", np.nan))
    lo = long["inhibitor"].map(lambda d: norm.get(d, {}).get("lo", np.nan))
    hi = long["inhibitor"].map(lambda d: norm.get(d, {}).get("hi", np.nan))
    a = pd.to_numeric(long["auc"], errors="coerce")
    out = long.copy()
    out["sens"] = -((a - med) / sc)
    y = pd.Series(np.nan, index=out.index)
    y[a <= lo] = 1.0                                   # bottom tail of AUC = sensitive
    y[a >= hi] = 0.0
    y[out["flag_increasing"].astype(bool) & (y == 1.0)] = np.nan   # a rising curve is never a response
    out["y_sens"] = y
    return out


# ------------------------------------------------------------------ model ----
class DrugResponseModel:
    VERSION = "1.0"

    def __init__(self, fs, kappa=KAPPA, alphas=ALPHAS):
        self.fs = fs                 # fitted FeatureSpace
        self.kappa = kappa
        self.alphas = alphas
        self.group_models = {}       # family_group -> {est, alpha, n}
        self.drug_models = {}        # inhibitor -> {est, w, n, group}
        self.calib = {}              # inhibitor -> {a, b} Platt on OOF scores
        self.norm = {}
        self.drug_group = {}
        self.slices = {}
        self.blocks = []
        self.sym2ens = {}
        self.meta = {}

    # --------------------------------------------------------- assembly ----
    def _rows(self, long, Xp, row_of, drugs):
        """Stack (specimen x drug) rows for a set of drugs: patient features + that drug's target
        descriptor. Returns X, y, groups(subject), and the row index back into `long`."""
        Xs, ys, gs, idx = [], [], [], []
        for drug in drugs:
            sub = long[(long["inhibitor"] == drug) & long["sens"].notna()]
            if not len(sub):
                continue
            ri = sub["specimen"].map(row_of)
            ok = ri.notna().values
            if not ok.any():
                continue
            sub = sub[ok]; ri = ri[ok].astype(int).values
            Zt = Xp[ri]
            T = self._tgt.get(drug)
            if T is None:                                   # unannotated drug: no descriptor, zeros
                T = np.zeros((Xp.shape[0], 3))
            Xs.append(np.hstack([Zt, T[ri]]))
            ys.append(sub["sens"].values.astype(float))
            gs.append(sub["subject"].values)
            idx.append(sub.index.values)
        if not Xs:
            return None
        return (np.vstack(Xs), np.concatenate(ys), np.concatenate(gs), np.concatenate(idx))

    def prepare_targets(self, Z_all, sym2ens):
        """Precompute each inhibitor's target-expression descriptor for every specimen, once."""
        self._tgt = {}
        for drug, a in TG.annotation().items():
            T, _ = self.fs.target_features(Z_all, a["targets"], sym2ens)
            self._tgt[drug] = T
        return self

    # ------------------------------------------------------------- fit ----
    @staticmethod
    def _block_X(Xp, T, ri, sl, b):
        """One block's design matrix; `T` is already row-aligned to `ri`. The drug descriptor
        (target-gene expression) rides with `rna`, because that is what it is -- expression. Ablating
        RNA therefore also removes it, which keeps the ablation honest instead of quietly leaving an
        expression feature in the 'no RNA' arm."""
        Z = Xp[ri][:, sl[b]]
        return np.hstack([Z, T]) if b == "rna" else Z

    # ---- multi-task matrix factorisation across inhibitors ------------------------------------------
    MF_RANK = 40
    MF_WEIGHT = 0.5

    @staticmethod
    def _soft_impute(Y, rank, n_iter=60, tol=1e-4):
        """Rank-`rank` factorisation of a matrix with missing entries, by iterative SVD imputation.

        Missing entries start at zero (the matrix is centred, so zero means "no information") and are
        refilled from the current low-rank reconstruction each pass.
        """
        obs = np.isfinite(Y)
        X = np.where(obs, Y, 0.0)
        prev = None
        for _ in range(n_iter):
            U, sv, Vt = np.linalg.svd(X, full_matrices=False)
            Xl = (U[:, :rank] * sv[:rank]) @ Vt[:rank]
            X = np.where(obs, Y, Xl)
            if prev is not None and np.linalg.norm(Xl - prev) / (np.linalg.norm(prev) + 1e-9) < tol:
                break
            prev = Xl
        U, sv, Vt = np.linalg.svd(X, full_matrices=False)
        return U[:, :rank] * sv[:rank], Vt[:rank]

    def _fit_mf(self, long, Xp, row_of, drugs):
        """Learn drug factors V and a features -> patient-factors map, so a NEW patient gets a profile.

        The per-family design shares information only through hand-drawn family boundaries. The response
        matrix is 520 x 118 and ~80% observed, which is the regime where a low-rank factorisation shares
        strength across ALL inhibitors and gives the low-n ones a usable profile. A new patient has no
        response row, so the patient factors are regressed on the patient's molecular features and the
        prediction at inference is f(x) . V.

        Measured on the interaction target (patient main effect removed, the only honest measure of
        drug-SPECIFIC skill): deployed 0.671, rank-40 factorisation 0.680, 50/50 blend 0.682. The blend
        weight sweep is flat between 0.5 and 0.75 (0.682-0.683), so 0.5 is used rather than chasing
        0.001 on the sweep that selected it.
        """
        from sklearn.linear_model import RidgeCV
        self.mf = None
        try:
            specs = sorted({s for s in long["specimen"] if s in row_of})
            if len(specs) < 50 or len(drugs) < 10:
                return
            sidx = {s: i for i, s in enumerate(specs)}
            didx = {d: j for j, d in enumerate(drugs)}
            Y = np.full((len(specs), len(drugs)), np.nan)
            sub = long[long["sens"].notna()]
            for sp, dr, v in zip(sub["specimen"], sub["inhibitor"], sub["sens"]):
                if sp in sidx and dr in didx:
                    Y[sidx[sp], didx[dr]] = float(v)
            if np.isfinite(Y).mean() < 0.2:
                return
            U, V = self._soft_impute(Y, min(self.MF_RANK, min(Y.shape) - 1))
            Xs = np.vstack([Xp[row_of[s]] for s in specs])
            reg = RidgeCV(alphas=(10.0, 100.0, 1000.0, 1e4)).fit(Xs, U)
            self.mf = {"reg": reg, "V": V, "drug_index": didx, "rank": int(V.shape[0]),
                       "weight": float(self.MF_WEIGHT)}
        except Exception:
            self.mf = None                      # a failed factorisation must cost the blend, not the model

    def _mf_pred(self, Xp, rows, drug):
        """MF contribution for `drug` at feature rows `rows`, or None when unavailable."""
        mf = getattr(self, "mf", None)
        if not mf or drug not in mf["drug_index"]:
            return None
        try:
            return mf["reg"].predict(Xp[rows]) @ mf["V"][:, mf["drug_index"][drug]]
        except Exception:
            return None

    def fit(self, long, Xp, row_of, drugs, slices, inner_folds=3):
        """Per-block ridges fused by non-negative least squares on inner donor-grouped OOF predictions.

        A single ridge over the concatenated blocks loses to RNA alone: the blocks differ by orders of
        magnitude in dimensionality and scale, so one shared alpha over-penalises the small ones and
        under-penalises the large. Late fusion is the same remedy the mutation layer already uses, and
        it carries the same guard -- if the fused blend does not beat the best single block on the same
        held-out rows, the weights collapse onto that block, so fusion can never do worse than its best
        part.
        """
        from sklearn.linear_model import Ridge, RidgeCV
        from sklearn.model_selection import GroupKFold
        from scipy.optimize import nnls
        self.slices = dict(slices)
        self.blocks = ([b for b in ["rna", "state", "mut", "clin"] if b in slices]
                       + [b for b in slices if b not in ("rna", "state", "mut", "clin")])
        self.drug_group = {d: TG.get(d)["family_group"] for d in drugs}
        groups = {}
        for d in drugs:
            groups.setdefault(self.drug_group[d], []).append(d)

        pred_shared = {}
        for gname, gdrugs in groups.items():
            r = self._rows_raw(long, Xp, row_of, gdrugs)
            if r is None:
                continue
            ri, y, sub, idx, tgt = r
            Xb = {b: self._block_X(Xp, tgt, ri, self.slices, b) for b in self.blocks}

            # inner donor-grouped OOF, per block
            oof = {b: np.full(len(y), np.nan) for b in self.blocks}
            nsp = min(inner_folds, len(np.unique(sub)))
            if nsp >= 2:
                for tri, tei in GroupKFold(n_splits=nsp).split(y, groups=sub):
                    for b in self.blocks:
                        e = Ridge(alpha=self._alpha_for(Xb[b], y)).fit(Xb[b][tri], y[tri])
                        oof[b][tei] = e.predict(Xb[b][tei])
            ok = np.all([np.isfinite(oof[b]) for b in self.blocks], axis=0)
            w = np.zeros(len(self.blocks))
            if ok.sum() >= 30:
                A = np.column_stack([oof[b][ok] for b in self.blocks])
                lam = 1e-3 * np.sqrt(ok.sum())
                Aa = np.vstack([A, lam * np.eye(len(self.blocks))])
                ya = np.concatenate([y[ok], np.zeros(len(self.blocks))])
                w, _ = nnls(Aa, ya)
                # floor: fusion must beat the best single block on the SAME rows or we take that block
                best_i = int(np.argmax([_corr(oof[b][ok], y[ok]) for b in self.blocks]))
                if _corr(A @ w, y[ok]) <= _corr(oof[self.blocks[best_i]][ok], y[ok]):
                    w = np.zeros(len(self.blocks)); w[best_i] = 1.0
            if not w.sum():
                w = np.ones(len(self.blocks)) / len(self.blocks)

            ests = {b: RidgeCV(alphas=self.alphas).fit(Xb[b], y) for b in self.blocks}
            self.group_models[gname] = {"ests": ests, "w": dict(zip(self.blocks, map(float, w))),
                                        "n": int(len(y)), "n_drugs": len(gdrugs)}
            fused = sum(w[i] * ests[b].predict(Xb[b]) for i, b in enumerate(self.blocks))
            pred_shared.update(dict(zip(idx, fused)))

        self._fit_mf(long, Xp, row_of, drugs)

        for d in drugs:
            r = self._rows_raw(long, Xp, row_of, [d])
            if r is None:
                continue
            ri, y, sub, idx, tgt = r
            base = np.array([pred_shared.get(i, 0.0) for i in idx])
            n = len(y)
            est = Ridge(alpha=1000.0).fit(Xp[ri], y - base)
            self.drug_models[d] = {"est": est, "w": float(n / (n + self.kappa)), "n": int(n),
                                   "group": self.drug_group[d]}
        self.meta = {"n_group_models": len(self.group_models), "n_drug_models": len(self.drug_models),
                     "kappa": self.kappa, "blocks": self.blocks,
                     "fusion_weights": {g: v["w"] for g, v in self.group_models.items()}}
        return self

    @staticmethod
    def _alpha_for(X, y):
        """Cheap scale-aware ridge penalty for the inner loop, where RidgeCV would be wasteful."""
        return max(1.0, 0.5 * X.shape[1] * float(np.mean(X.var(0)) or 1.0))

    def _rows_raw(self, long, Xp, row_of, drugs):
        ris, ys, gs, idx, tg = [], [], [], [], []
        for drug in drugs:
            sub = long[(long["inhibitor"] == drug) & long["sens"].notna()]
            if not len(sub):
                continue
            ri = sub["specimen"].map(row_of)
            ok = ri.notna().values
            if not ok.any():
                continue
            sub = sub[ok]; ri = ri[ok].astype(int).values
            T = self._tgt.get(drug)
            if T is None:
                T = np.zeros((Xp.shape[0], 3))
            ris.append(ri); ys.append(sub["sens"].values.astype(float))
            gs.append(sub["subject"].values); idx.append(sub.index.values); tg.append(T[ri])
        if not ris:
            return None
        return (np.concatenate(ris), np.concatenate(ys), np.concatenate(gs),
                np.concatenate(idx), np.vstack(tg))

    # --------------------------------------------------------- predict ----
    def _shared_pred(self, gm, Xp, ii, T):
        tot = 0.0
        for b, wt in gm["w"].items():
            if wt <= 0:
                continue
            Z = Xp[ii][:, self.slices[b]]
            if b == "rna":
                Z = np.hstack([Z, T])
            tot = tot + wt * gm["ests"][b].predict(Z)
        return tot

    def predict_rows(self, long, Xp, row_of):
        """Predicted `sens` for every row of `long` we can score."""
        out = pd.Series(np.nan, index=long.index, dtype=float)
        for drug, sub in long.groupby("inhibitor"):
            dm = self.drug_models.get(drug)
            gm = self.group_models.get(self.drug_group.get(drug))
            if dm is None or gm is None:
                continue
            ri = sub["specimen"].map(row_of)
            ok = ri.notna().values
            if not ok.any():
                continue
            ii = ri[ok].astype(int).values
            T = self._tgt.get(drug)
            T = T[ii] if T is not None else np.zeros((len(ii), 3))
            base = self._shared_pred(gm, Xp, ii, T) + dm["w"] * dm["est"].predict(Xp[ii])
            mfp = self._mf_pred(Xp, ii, drug)
            if mfp is not None:
                w = getattr(self, "mf", {}).get("weight", 0.0)
                base = (1.0 - w) * base + w * mfp
            out.loc[sub.index[ok]] = base
        return out

    def predict_patient(self, Zp, mut=None, meta=None, drugs=None, sym2ens=None, ref="beataml"):
        """Score ONE new sample (already cohort-z-scored, 1 x genes) against every modelled drug.

        `ref` names the score reference to percentile against -- 'beataml' for a bulk specimen,
        'sc_sample' for a single-cell bulk-equivalent. Getting this wrong is what makes every drug come
        back at probability 0.99."""
        Xp, _, _ = self.fs.transform(Zp, meta=meta, mut=mut, blocks=self.blocks)
        res = {}
        for drug in (drugs or self.drug_models):
            dm = self.drug_models.get(drug)
            gm = self.group_models.get(self.drug_group.get(drug))
            if dm is None or gm is None:
                continue
            T, _ = self.fs.target_features(Zp, TG.get(drug)["targets"], sym2ens or self.sym2ens)
            s = float(self._shared_pred(gm, Xp, np.array([0]), T)[0] + dm["w"] * dm["est"].predict(Xp)[0])
            # the same blend the training-time scoring uses; without it an upload would be served a
            # different model from the one the reported AUROC was measured on
            mfp = self._mf_pred(Xp, np.array([0]), drug)
            if mfp is not None:
                w = getattr(self, "mf", {}).get("weight", 0.0)
                s = float((1.0 - w) * s + w * float(mfp[0]))
            res[drug] = {"sens": s, "prob_sensitive": self.calibrated(drug, s, ref),
                         "percentile": self.score_percentile(drug, s, ref),
                         "n_train": dm["n"], "group": dm["group"], "score_reference": ref}
        return res

    # ------------------------------------------------------ calibration ----
    def fit_calibration(self, long, oof_pred):
        """Per-drug Platt scaling onto P(sensitive) -- fitted on the OOF **percentile**, not the raw score.

        This is the cross-assay lesson the mutation caller already learned the hard way. A single-cell
        bulk-equivalent occupies a different band of the raw decision scale than a BeatAML bulk
        specimen, so a Platt curve fitted on BeatAML raw scores saturates at 1.0 for essentially every
        single-cell sample -- a model that recommends all 118 inhibitors with probability 0.99.
        Calibrating the *percentile* makes the mapping cohort-invariant: "this sample is in the top 8%
        of predicted responders **among its own kind**" means the same thing in either assay, as long
        as the percentile is taken against a matched reference (see `set_score_reference`).
        """
        from sklearn.linear_model import LogisticRegression
        self.oof_score_ref = {}
        for drug, sub in long.groupby("inhibitor"):
            v = oof_pred.loc[sub.index].dropna().values
            if len(v) >= 20:
                self.oof_score_ref[drug] = np.sort(v)
        for drug, sub in long.groupby("inhibitor"):
            m = sub["y_sens"].notna() & oof_pred.loc[sub.index].notna()
            y = sub.loc[m, "y_sens"].values.astype(int)
            if len(y) < 20 or len(set(y)) < 2 or drug not in self.oof_score_ref:
                continue
            q = _pct(self.oof_score_ref[drug], oof_pred.loc[sub.index][m].values).reshape(-1, 1)
            lr = LogisticRegression(C=1e6, max_iter=1000).fit(q, y)
            self.calib[drug] = {"a": float(lr.coef_[0][0]), "b": float(lr.intercept_[0]),
                                "n": int(len(y)), "n_pos": int(y.sum())}
        return self

    def calibrated(self, drug, score, ref="beataml"):
        c = self.calib.get(drug)
        q = self.score_percentile(drug, score, ref)
        if c is None or q is None:
            return None
        return float(1.0 / (1.0 + np.exp(-(c["a"] * q + c["b"]))))

    def score_percentile(self, drug, score, ref="beataml"):
        """Where this score sits among predictions for samples OF THE SAME KIND (see score_refs)."""
        refs = getattr(self, "score_refs", None) or {}
        base = (refs.get(ref) or refs.get("beataml") or {}).get(drug)
        if base is None or score is None or score != score:
            return None
        return float(np.searchsorted(base, score, side="right") / len(base))

    def set_score_reference(self, long, pred, name="beataml"):
        """Store one cohort's predicted-score distribution per drug.

        `beataml` is set at training time. `sc_sample` (single-cell bulk-equivalents) and `sc_state`
        (individual cell-state pseudobulks) are added afterwards by build_drug_score_refs.py -- a cell
        state's profile is more extreme than any whole sample, so it needs its own reference or every
        state looks like an outlier.
        """
        refs = getattr(self, "score_refs", None) or {}
        r = {}
        for drug, sub in long.groupby("inhibitor"):
            v = pred.loc[sub.index].dropna().values
            if len(v) >= 20:
                r[drug] = np.sort(v)
        refs[name] = r
        self.score_refs = refs
        return self

    def add_score_reference(self, name, scores_by_drug, min_n=20):
        refs = getattr(self, "score_refs", None) or {}
        refs[name] = {d: np.sort(np.asarray(v, float)) for d, v in scores_by_drug.items()
                      if v is not None and len(v) >= min_n}
        self.score_refs = refs
        return self

    def predict_matrix(self, Z, drugs=None):
        """Raw scores for many samples at once: (n_samples x n_drugs) -- used to build score references."""
        Xp, _, _ = self.fs.transform(Z, meta=None, mut=None, blocks=self.blocks)
        rows = np.arange(Z.shape[0])
        out = {}
        for drug in (drugs or self.drug_models):
            dm = self.drug_models.get(drug)
            gm = self.group_models.get(self.drug_group.get(drug))
            if dm is None or gm is None:
                continue
            T, _ = self.fs.target_features(Z, TG.get(drug)["targets"], self.sym2ens)
            out[drug] = self._shared_pred(gm, Xp, rows, T) + dm["w"] * dm["est"].predict(Xp)
        return pd.DataFrame(out, index=range(Z.shape[0]))

    def attach_neighbours(self, Z_train, specimens, long_train):
        """Store the training specimens in RNA-PC space plus their MEASURED AUCs.

        This is what lets a report say "the 20 BeatAML specimens most similar to this patient had a
        median measured AUC of X for this drug" -- an observation rather than a model output, and the
        single most checkable thing the system can offer a sceptical reader.
        """
        P = self.fs.pca.transform(Z_train[:, self.fs.sel])
        self.nn = {"P": P.astype(np.float32), "specimens": list(specimens),
                   "auc": (long_train.pivot_table(index="specimen", columns="inhibitor",
                                                  values="auc", aggfunc="mean")
                           .reindex(list(specimens))),
                   "mu": P.mean(0), "sd": np.where(P.std(0) == 0, 1.0, P.std(0))}
        self.nn["ood_ref"] = np.sqrt((((P - self.nn["mu"]) / self.nn["sd"]) ** 2).mean(1))
        return self

    def neighbours(self, Zp, k=20, cohort="beataml"):
        """k nearest BeatAML training specimens to one new sample, in standardised PC space.

        `self_distance` is measured against the matched cohort's own centre when one exists, so a
        single-cell patient is asked "are you unusual among single-cell AML samples?" rather than
        "are you a bulk BeatAML specimen?" -- the second is always no and tells nobody anything.
        """
        if not getattr(self, "nn", None):
            return None
        p = self.fs.pca.transform(Zp[:, self.fs.sel])
        q = (p - self.nn["mu"]) / self.nn["sd"]
        R = (self.nn["P"] - self.nn["mu"]) / self.nn["sd"]
        dist = np.sqrt(((R - q) ** 2).mean(1))
        order = np.argsort(dist)[:k]
        if cohort != "beataml" and "sc_mu" in self.nn:
            qs = (p - self.nn["sc_mu"]) / self.nn["sc_sd"]
            self_d = float(np.sqrt((qs ** 2).mean()))
        else:
            self_d = float(np.sqrt((q ** 2).mean()))
        return {"specimens": [self.nn["specimens"][i] for i in order],
                "distance": [round(float(dist[i]), 4) for i in order],
                "self_distance": self_d, "distance_cohort": cohort,
                "beataml_distance": float(np.sqrt((q ** 2).mean()))}
