"""Deployable PRIMARY mutation caller — variant-level, from BULK RNA (or a single-cell bulk-equivalent).

Trained on BeatAML2 (n=707, full WES) over ~50 variant-level categories (FLT3_ITD / TKD, DNMT3A_R882 vs
nonR882, NRAS/KRAS G12/G13/Q61, TP53 hotspot-DBD vs LOF, splicing factors, ...) — the coverage the
single-cell multimodal system can't reach. One clean linear model per category on the 14,237 genes shared
across BeatAML / Leucegene / our atlas, so the SAME model scores plain bulk RNA-seq AND a single-cell
sample collapsed to its bulk-equivalent.

Cross-cohort calibration: every cohort is put on a common log2 scale then z-scored per cohort (handles the
platform/normalization batch difference). At predict time a new sample is z-scored against a stored
reference distribution (`sc` for atlas bulk-equivalents, `beataml` for bulk RNA-seq), then scored by the
BeatAML-trained model and percentile-mapped to [0,1] against the training score distribution.

  Train:   BulkMutationPredictor.train_from_bundle("bundle_data.npz", sym2ens) -> bulk_mutation_predictor.pkl
  Predict: P.predict(bulk_series, ref="sc")  # series indexed by ENSG or gene symbol -> {cat: {prob, call, ...}}
"""
from __future__ import annotations
import pickle
import numpy as np


def _pct(sorted_train, v):
    n = len(sorted_train); v = np.asarray(v, float); out = np.full(len(v), np.nan); ok = ~np.isnan(v)
    if n >= 2:
        out[ok] = np.searchsorted(sorted_train, v[ok], side="right") / n
    elif n == 1:
        out[ok] = 0.5
    return out


class BulkMutationPredictor:
    VERSION = "1.0"

    def __init__(self):
        self.genes = []            # ENSG order of the feature space (14,237)
        self.refs = {}             # ref_name -> (mu[genes], sd[genes]) on the common log2 scale
        self.score_refs = {}       # ref_name -> {cat: sorted decision scores OVER THAT COHORT}
        self.categories = []       # variant-level category names
        self.models = {}           # cat -> dict(est, sel, sorted_scores, sign, threshold, cv_auroc, n_pos)
        self.sym2ens = {}          # gene symbol -> ENSG (to accept symbol-keyed samples)
        self.meta = {}

    # ---------- helpers ----------
    @staticmethod
    def _clog(x_lin):
        return np.log2(np.clip(np.asarray(x_lin, float), 0, None) + 1.0)

    def _align(self, sample):
        """sample: dict / pandas Series keyed by ENSG or gene symbol -> linear vector over self.genes."""
        if not hasattr(self, "_gidx"):
            self._gidx = set(self.genes)
        if hasattr(sample, "to_dict"):
            sample = sample.to_dict()
        g2v = {}
        for k, v in sample.items():
            k = str(k)
            if k in self._gidx:
                g2v[k] = v
            else:
                e = self.sym2ens.get(k)
                if e in self._gidx:
                    g2v[e] = v
        return np.array([float(g2v.get(g, 0.0)) for g in self.genes], float)

    def _z(self, clog, ref):
        mu, sd = self.refs[ref]
        return (clog - mu) / sd

    # ---------- inference ----------
    def predict_one(self, cat, z, ref="sc"):
        """Percentile the sample against ITS OWN cohort's score distribution.

        Using the BeatAML score distribution for a non-BeatAML sample is wrong: sc scores occupy a
        compressed band of that scale, so high thresholds become literally unreachable (SRSF2 ceilinged
        at 0.883 vs a 0.888 cut -> 0/387 sc samples callable). With a cohort-matched reference the
        percentile means "where this sample sits among its own cohort", and a threshold means
        "call the top (1-thr) fraction" in any cohort.
        """
        m = self.models[cat]
        d = m["est"].decision_function(z[m["sel"]].reshape(1, -1))
        base = (self.score_refs.get(ref) or {}).get(cat)
        if base is None:
            base = m["sorted_scores"]                 # back-compat: pre-score_refs models
        q = _pct(base, d)[0]
        if q != q:
            return None
        prob = q if m["sign"] > 0 else 1.0 - q
        thr = m["threshold"]; ca = m["cv_auroc"]
        conf = "ok" if (ca is not None and ca >= 0.65) else ("abstain: weak (CV AUROC %.2f)" % ca if ca is not None else "abstain: unknown")
        return {"probability": round(float(prob), 3), "call": "present" if prob >= thr else "absent",
                "threshold": round(float(thr), 3), "cv_auroc": ca, "n_pos_train": m["n_pos"], "confidence": conf}

    def predict(self, sample, ref="sc"):
        if not hasattr(self, "_gidx"):
            self._gidx = set(self.genes)
        if ref not in self.refs:
            ref = next(iter(self.refs))
        z = self._z(self._clog(self._align(sample)), ref)
        out = {}
        for cat in self.categories:
            r = self.predict_one(cat, z, ref)         # percentile against the SAME cohort we z-scored to
            if r is not None:
                out[cat] = r
        return dict(sorted(out.items(), key=lambda kv: -(kv[1]["probability"] or 0)))

    # ---------- persistence ----------
    def save(self, path):
        with open(path, "wb") as f:
            pickle.dump(self, f)

    @staticmethod
    def load(path):
        with open(path, "rb") as f:
            return pickle.load(f)

    def summary(self):
        aucs = [m["cv_auroc"] for m in self.models.values() if m["cv_auroc"] is not None]
        return {"version": self.VERSION, "n_categories": len(self.categories),
                "mean_cv_auroc": round(float(np.mean(aucs)), 3) if aucs else None,
                "refs": list(self.refs), "genes": len(self.genes), **self.meta}

    # ---------- training ----------
    @classmethod
    def train_from_bundle(cls, bundle_path, sym2ens, varcap=2500, min_pos=6, base="logL2"):
        from sklearn.linear_model import LogisticRegression
        from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
        from sklearn.model_selection import StratifiedKFold
        from sklearn.metrics import roc_auc_score
        d = np.load(bundle_path, allow_pickle=True)
        genes = [str(g) for g in d["genes"]]; cats = [str(c) for c in d["drivers"]]
        baX = d["ba_X"].astype(float); baL = d["ba_L"].astype(float)

        def clog_z_ref(X):
            cl = np.log2(np.clip(X, 0, None) + 1.0)
            mu = cl.mean(0); sd = cl.std(0); sd[sd == 0] = 1.0
            return mu, sd, (cl - mu) / sd

        self = cls(); self.genes = genes; self.sym2ens = dict(sym2ens)
        Zs = {}
        for name, key in [("beataml", "ba_X"), ("sc", "sc_X"), ("leucegene", "lg_X")]:
            mu, sd, Z = clog_z_ref(d[key].astype(float))
            self.refs[name] = (mu, sd); Zs[name] = Z          # keep each cohort's z-matrix for score refs
        Zba = Zs["beataml"]

        def mk(name):
            if name == "shrLDA":
                return LinearDiscriminantAnalysis(solver="lsqr", shrinkage=0.5)
            return LogisticRegression(C=1.0, class_weight="balanced", max_iter=2000)

        def topvar(Z):
            return np.argsort(Z.var(0))[::-1][:varcap]

        for j, cat in enumerate(cats):
            y = baL[:, j]; ok = ~np.isnan(y); yv = y[ok].astype(int); Z = Zba[ok]
            if yv.sum() < min_pos or (yv == 0).sum() < min_pos:
                continue
            # 5-fold CV-OOF for AUROC + F1-max threshold + orientation
            oof = np.full(len(yv), np.nan)
            ns = min(5, int(yv.sum()), int((yv == 0).sum()))
            for tri, tei in StratifiedKFold(ns, shuffle=True, random_state=0).split(Z, yv):
                sel = topvar(Z[tri]); est = mk(base).fit(Z[tri][:, sel], yv[tri])
                oof[tei] = est.decision_function(Z[tei][:, sel])
            om = ~np.isnan(oof)
            try:
                a = roc_auc_score(yv[om], oof[om])
            except Exception:
                a = np.nan
            sign = 1.0 if (a != a or a >= 0.5) else -1.0
            cv_auroc = round(float(max(a, 1 - a)), 3) if a == a else None
            # deployed model: refit on ALL BeatAML
            sel = topvar(Z); est = mk(base).fit(Z[:, sel], yv)
            dtr = est.decision_function(Z[:, sel]); ssorted = np.sort(dtr)
            # F1-max threshold on the OOF percentiles (oriented)
            p = _pct(np.sort(oof[om]), oof[om]); p = p if sign > 0 else 1 - p
            yv_oof = yv[om]; bt, bf1 = 0.5, -1.0
            for t in np.unique(p):
                pr = p >= t
                tp = int((pr & (yv_oof == 1)).sum()); fp = int((pr & (yv_oof == 0)).sum()); fn = int(((~pr) & (yv_oof == 1)).sum())
                prec = tp / (tp + fp) if (tp + fp) else 0.0; rec = tp / (tp + fn) if (tp + fn) else 0.0
                f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
                if f1 > bf1:
                    bf1, bt = f1, float(t)
            # cohort-matched score references: each cohort's own decision-score distribution, so a
            # sample is percentiled among its peers rather than against BeatAML's scale.
            for _name, _Z in Zs.items():
                self.score_refs.setdefault(_name, {})[cat] = np.sort(est.decision_function(_Z[:, sel]))
            self.categories.append(cat)
            self.models[cat] = {"est": est, "sel": sel, "sorted_scores": ssorted, "sign": sign,
                                "threshold": round(bt, 3), "cv_auroc": cv_auroc, "n_pos": int(yv.sum())}
        self.meta = {"trained_on": "BeatAML2 (n=707, WES)", "base_model": base, "varcap": varcap,
                     "min_pos": min_pos}
        return self
