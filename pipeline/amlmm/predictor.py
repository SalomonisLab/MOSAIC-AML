"""Deployable per-mutation MOSAIC-AML predictor — the validated combiner, consolidated.

Per mutation, each modality has a frozen pipeline: StandardScaler -> differential feature selection ->
LinearSVC. Its decision scores are percentile-normalised against the TRAIN distribution and oriented;
the per-modality scores are then combined with per-mutation optimised modality weights (ridge-NNLS,
learned on donor-grouped OOF — the "molded weights" approach, sealed held-out mean AUC ~0.86, and
provably >= the best single modality on the data it is fit on).

Late fusion by design: at predict time the weights renormalise over whatever modalities a sample
actually has, so an upload with only Composition + RNA degrades gracefully instead of breaking.

  Train once on the cohort (see pipeline/train_predictor.py) -> mutation_predictor.pkl
  Load + predict anywhere (sklearn + numpy only):
      P = MutationPredictor.load("mutation_predictor.pkl")
      calls = P.predict({"RNA": series_or_vec, "Composition": ...})   # -> {mut: {prob, call, contributions}}
"""
from __future__ import annotations
import pickle
import numpy as np


def diff_select(X, y, k):
    """top-(k/2) up + top-(k/2) down differential features (train-only); leakage-safe per call."""
    from sklearn.feature_selection import f_classif
    if X.shape[1] <= k:
        return np.arange(X.shape[1])
    F = np.nan_to_num(f_classif(X, y)[0]); md = X[y == 1].mean(0) - X[y == 0].mean(0)
    o = np.argsort(np.sign(md) * F)
    return np.unique(np.concatenate([o[:k // 2], o[-(k // 2):]]))


def _pct(sorted_train, v):
    n = len(sorted_train)
    v = np.asarray(v, float); out = np.full(len(v), np.nan); ok = ~np.isnan(v)
    if n >= 2:
        out[ok] = np.searchsorted(sorted_train, v[ok], side="right") / n
    elif n == 1:
        out[ok] = 0.5
    return out


class ModalityModel:
    """Frozen per-(mutation, modality) pipeline + the train score distribution for percentile mapping."""
    __slots__ = ("features", "scaler", "keep", "sel", "svm", "sorted_scores", "sign", "oof_auc")

    def __init__(self, features, scaler, keep, sel, svm, sorted_scores, sign, oof_auc):
        self.features = list(features)          # modality's full feature names (to align a new sample)
        self.scaler = scaler; self.keep = keep; self.sel = sel; self.svm = svm
        self.sorted_scores = sorted_scores; self.sign = float(sign); self.oof_auc = oof_auc

    def _align(self, sample):
        """sample: 1-D vector aligned to self.features, OR a {feature_name: value} mapping. Mapping lookup
        is prefix-robust: a 'comp::ASDC' feature matches a plain 'ASDC' key and vice-versa, so callers
        (validation, ingest_patient) can pass raw cell-state / gene names without knowing the prefix."""
        if hasattr(sample, "get") and not isinstance(sample, (list, tuple, np.ndarray)):
            out = []
            for f in self.features:
                v = sample.get(f)
                if v is None:
                    v = sample.get(f.split("::", 1)[-1])           # feature is prefixed, key is not
                if v is None and "::" not in str(f):
                    v = sample.get("comp::" + f)                   # key is prefixed, feature is not
                out.append(float(v) if v is not None else 0.0)
            return np.array(out, float)
        v = np.asarray(sample, float)
        return v if v.shape[0] == len(self.features) else np.resize(v, len(self.features))

    def score(self, sample):
        """Return the oriented percentile in [0,1] (higher = more likely positive)."""
        x = self._align(sample).reshape(1, -1)
        xk = self.scaler.transform(x[:, self.keep])[:, self.sel]
        d = self.svm.decision_function(xk)
        q = _pct(self.sorted_scores, d)[0]
        if q != q:
            return np.nan
        return q if self.sign > 0 else 1.0 - q


class MutationPredictor:
    VERSION = "1.0"

    def __init__(self):
        self.mutations = []          # mutation flags (mut_* / cyto_*)
        self.modalities = []         # modality names in train order
        self.models = {}             # (mutation, modality) -> ModalityModel
        self.weights = {}            # mutation -> {modality: weight}  (sum 1, deployable)
        self.heldout_auc = {}        # mutation -> sealed held-out AUC (display)
        self.train_auc = {}          # mutation -> donor-grouped CV (out-of-fold) AUROC on the TRAINING set
        self.prevalence = {}         # mutation -> cohort positive rate
        self.thresholds = {}         # mutation -> calibrated present/absent cut (cohort Youden-J; default 0.5)
        self.meta = {}               # rep notes (rna='raw+prog', n_train, date, ...)

    # ---- persistence ----
    def save(self, path):
        with open(path, "wb") as f:
            pickle.dump(self, f)

    @staticmethod
    def load(path):
        with open(path, "rb") as f:
            return pickle.load(f)

    # ---- inference ----
    def predict_one(self, mutation, sample):
        """sample: {modality: vector|mapping}. Returns {prob, call, contributions, n_modalities}."""
        w = self.weights.get(mutation, {})
        contribs, num, den = [], 0.0, 0.0
        for mod in self.modalities:
            mm = self.models.get((mutation, mod))
            if mm is None or mod not in sample or sample[mod] is None:
                continue
            s = mm.score(sample[mod])
            if s != s:
                continue
            wm = float(w.get(mod, 0.0))
            contribs.append({"modality": mod, "weight": round(wm, 3), "score": round(float(s), 3),
                             "oof_auc": (round(mm.oof_auc, 3) if mm.oof_auc is not None else None)})
            if wm > 0:
                num += wm * s; den += wm
        prob = (num / den) if den > 0 else (np.mean([c["score"] for c in contribs]) if contribs else None)
        contribs.sort(key=lambda c: -c["weight"])
        thr = float(getattr(self, "thresholds", {}).get(mutation, 0.5))     # calibrated cut (over-calling fix)
        return {"probability": (round(float(prob), 3) if prob is not None else None),
                "call": (None if prob is None else ("present" if prob >= thr else "absent")),
                "threshold": round(thr, 3),
                "n_modalities": len(contribs), "contributions": contribs,
                "heldout_auc": self.heldout_auc.get(mutation)}

    def predict(self, sample, mutations=None):
        """Per-mutation predictions for one sample. -> {mutation: predict_one(...)}, sorted by probability."""
        muts = mutations or self.mutations
        out = {m: self.predict_one(m, sample) for m in muts}
        return dict(sorted(out.items(), key=lambda kv: -(kv[1]["probability"] or 0)))

    def summary(self):
        return {"version": self.VERSION, "n_mutations": len(self.mutations),
                "modalities": self.modalities, "mean_heldout_auc":
                    (round(float(np.nanmean([v for v in self.heldout_auc.values() if v is not None])), 3)
                     if self.heldout_auc else None), **self.meta}
