"""Cohort adapter — route each mutation category to whichever model is better for ONE cohort.

The shipped BulkMutationPredictor is trained on BeatAML (n=707). On a cohort from another platform it
can rank some drivers well and others at chance (dataset 32266596, BD Rhapsody: NPM1 0.94 but RUNX1
0.53). A model fitted on that cohort fixes some drivers and breaks others: common drivers gain, rare
ones (IDH1, IDH2, TET2, FLT3-TKD at 9-16 positives) lose to the 707-patient model. This module keeps
both and lets each category use the one that is better *for this cohort*:

  * a unit (a driver) is ROUTED to the cohort model only when, inside training folds, the cohort model
    beats the shipped model by a margin on ranking (AUROC) AND is not worse at its operating point
    (F1). Otherwise it stays on the shipped model, untouched. The choice is made with training data
    only; `train_cohort_adapter.py` measures the whole procedure by nested CV and on a sealed split.
  * two cohort-model arms compete per unit: `reweight` (a logistic layer over the shipped model's own
    50 category scores - a warm start that keeps BeatAML's knowledge) and `scratch` (the bulk caller's
    recipe refitted on this cohort's expression).
  * units the shipped caller has no category for (cytogenetics) are ADDED only if their cross-validated
    AUROC clears a floor.

It is opt-in and cohort-specific. `apply` refuses (returns the shipped calls unchanged, with a note)
when the input does not cover the adapter's genes or does not look like the cohort it was fitted on,
so attaching the wrong adapter cannot silently rewrite calls for another platform.

The adapter file is plain JSON (no pickle): coefficients, per-gene reference statistics, and each
unit's out-of-fold score distribution. It holds no per-patient expression.
"""
from __future__ import annotations
import gzip
import json
import numpy as np

VERSION = "1.0"          # format written by fit() below; v2 files are written by the custom-model export
LN2 = np.log(2.0)


def _pct(sorted_ref, v):
    ref = np.asarray(sorted_ref, float)
    if ref.size == 0:
        return float("nan")
    return float(np.searchsorted(ref, v, side="right") / ref.size)


def _cp10k(x):
    x = np.clip(np.asarray(x, float), 0, None)
    s = x.sum()
    return x / s * 1e4 if s > 0 else x


def shipped_decisions(bp, sample, ref="sc", mu_sd=None):
    """The shipped model's oriented decision score for every category.

    ref: one of the shipped model's own z-references. mu_sd: an explicit per-gene (mean, sd) in the model's
    log2 space instead - a COHORT-MATCHED reference. On dataset 32266596 that label-free switch alone took
    the shipped model's mean AUROC from 0.68 to 0.80 on held-out patients."""
    clog = bp._clog(bp._align(sample))
    z = bp._z(clog, ref) if mu_sd is None else (clog - mu_sd[0]) / mu_sd[1]
    out = np.empty(len(bp.categories))
    for j, c in enumerate(bp.categories):
        m = bp.models[c]
        out[j] = float(m["est"].decision_function(z[m["sel"]].reshape(1, -1))[0]) * m["sign"]
    return out


class CohortAdapter:
    """Two file formats share this class.

    v1 (`train_cohort_adapter.py`): per unit one of two arms, `reweight` (a layer over the shipped model's
        50 scores) or `scratch` (bulk-RNA genes).
    v2 (custom per-mutation models): per unit ONE linear model over named input blocks. A unit's ensemble
        of linear candidates is collapsed into `coef[block]` + intercept, so scoring is
            score = sum_blocks z_block . coef_block + intercept ,   z = (feature - mean) / sd.
        Block kinds: `genes` (log2 CP10k), `shipped_scores` (the shipped model's 50 decision scores) and
        `dosage` (expression dosage per chromosome arm, relative to the cohort's per-gene median).
        A v2 file may also carry `cohort_reference`: per-gene mean/sd of the cohort in the shipped model's
        space plus each category's score distribution. Every shipped category is then re-referenced to the
        cohort first (label-free), and the custom models are written on top of that.
    """

    def __init__(self, d):
        self.d = d
        self.cohort = d["cohort"]
        self.v2 = str(d.get("version", "1.0")).startswith("2")
        g = d["guard_profile"] if self.v2 else d["features"]
        self.genes = list(g["genes"])
        self.g_mu = np.asarray(g["mu"], float)
        self.g_sd = np.asarray(g["sd"], float)
        if not self.v2:
            self.dec_mu = np.asarray(d["features"]["dec_mu"], float)
            self.dec_sd = np.asarray(d["features"]["dec_sd"], float)
        self.units = d["units"]
        cr = d.get("cohort_reference")
        self.cref = None if not cr else (np.asarray(cr["mu"], float), np.asarray(cr["sd"], float))

    # ------------------------------------------------------------------ io ----
    @classmethod
    def load(cls, path):
        op = gzip.open if str(path).endswith(".gz") else open
        with op(path, "rt", encoding="utf-8") as f:
            return cls(json.load(f))

    def save(self, path):
        op = gzip.open if str(path).endswith(".gz") else open
        with op(path, "wt", encoding="utf-8") as f:
            json.dump(self.d, f, indent=None, separators=(",", ":"))

    def summary(self):
        routed = [u for u in self.units if u["routed"]]
        return {"cohort": self.cohort, "version": self.d.get("version"),
                "routed": [{"unit": u["name"], "arm": u.get("arm", "custom"), "mode": u["mode"],
                            "cohort_cv_auroc": u["cv_auroc"], "shipped_auroc": u.get("shipped_auroc")}
                           for u in routed],
                "kept_on_shipped": [{"unit": u["name"], "why": u["reason"]} for u in self.units if not u["routed"]],
                "evaluation": self.d.get("evaluation"),
                "note": ("categories listed under `routed` are called by a model fitted on cohort %s; every "
                         "other category is the shipped BeatAML-trained call, unchanged" % self.cohort)}

    # ------------------------------------------------------------ features ----
    @staticmethod
    def _linear(bp, sample):
        """gene -> CP10k lookup for one sample, renormalised over everything it carries (so counts, CP10k
        and TPM inputs land on the same scale), with a symbol->ENSG fallback."""
        s = sample.to_dict() if hasattr(sample, "to_dict") else dict(sample)
        keys = list(s)
        lin = dict(zip(keys, _cp10k([float(s[k]) for k in keys])))
        s2e = getattr(bp, "sym2ens", {}) or {}

        def get(g):
            v = lin.get(g)
            if v is None:
                e = s2e.get(g)
                v = lin.get(e) if e else None
            return v
        return get

    @staticmethod
    def _log_vector(get, genes):
        x = np.full(len(genes), np.nan)
        for i, g in enumerate(genes):
            v = get(g)
            if v is not None:
                x[i] = np.log1p(v) / LN2
        return x

    def _gene_vector(self, bp, sample):
        """log2(CP10k+1) over the guard genes. Returns (x, coverage)."""
        x = self._log_vector(self._linear(bp, sample), self.genes)
        return x, float(np.isfinite(x).mean())

    def reference_fit(self, x):
        z = (x - self.g_mu) / self.g_sd
        z = z[np.isfinite(z)]
        return float(np.mean(np.abs(z))) if z.size else float("inf")

    def profile_corr(self, x):
        """Correlation of the sample's expression profile with the cohort's mean profile. This is the
        guard that separates platforms: mean|z| does not (32266596 reaches 1.30, 10x atlas samples start
        at 1.16), but the profile correlation does (32266596 median 0.93; 10x samples <= 0.84)."""
        return float(np.corrcoef(x, self.g_mu)[0, 1])

    def _decisions(self, bp, sample, reference):
        """Shipped decision scores under the shipped reference or the cohort-matched one, computed once per sample."""
        if reference not in self._dec_cache:          # cache is reset at the start of every apply()/score_units()
            mu_sd = self.cref if reference == "cohort" else None
            self._dec_cache[reference] = shipped_decisions(bp, sample, self.d["shipped_ref"], mu_sd)
        return self._dec_cache[reference]

    def _block_z(self, name, bp, sample, get):
        """Standardised features of one v2 block for one sample (missing -> 0, i.e. the cohort mean)."""
        b = self.d["blocks"][name]
        if b["kind"] == "genes":
            x = self._log_vector(get, b["genes"])
        elif b["kind"] == "shipped_scores":
            x = self._decisions(bp, sample, b.get("reference", "shipped"))
        elif b["kind"] == "dosage":
            x = dosage_features(self._log_vector(get, b["definition"]["genes"]), b["definition"])
        else:
            raise ValueError("unknown block kind %r" % b["kind"])
        if "sel" in b:                                   # features kept when the block was fitted
            x = x[np.asarray(b["sel"], int)]
        return np.nan_to_num((x - np.asarray(b["mu"], float)) / np.asarray(b["sd"], float))

    # --------------------------------------------------------------- apply ----
    def apply(self, bp, sample, calls):
        """calls: the shipped {category: call-dict}. Returns a new dict with routed units replaced and
        added units appended. Never mutates the input."""
        out = {k: dict(v) for k, v in calls.items()}
        self._dec_cache = {}
        get = self._linear(bp, sample)
        x = self._log_vector(get, self.genes)
        cov = float(np.isfinite(x).mean())
        lim = self.d["guards"]
        xf = np.where(np.isfinite(x), x, self.g_mu)
        mz = self.reference_fit(xf)
        pc = self.profile_corr(xf)
        refuse = None
        if cov < lim["min_gene_coverage"]:
            refuse = "only %.0f%% of the adapter's genes are present in this sample (need %.0f%%)" % (
                100 * cov, 100 * lim["min_gene_coverage"])
        elif pc < lim.get("min_profile_corr", -1.0):
            refuse = ("expression profile correlates %.3f with cohort %s (limit %.3f): most likely a different "
                      "platform or tissue than the adapter was fitted on" % (pc, self.cohort, lim["min_profile_corr"]))
        elif mz > lim["max_mean_abs_z"]:
            refuse = ("mean|z| %.2f against cohort %s (limit %.2f): this sample does not look like the cohort "
                      "the adapter was fitted on" % (mz, self.cohort, lim["max_mean_abs_z"]))
        if refuse:
            for r in out.values():
                r.setdefault("cohort_adapter", {"applied": False, "cohort": self.cohort, "reason": refuse})
            return out
        prov = {"applied": True, "cohort": self.cohort, "mean_abs_z": round(mz, 3), "profile_corr": round(pc, 3)}
        if self.v2:
            if self.cref is not None:            # layer 1: every shipped category, re-referenced to the cohort
                dec = self._decisions(bp, sample, "cohort")
                sref = self.d["cohort_reference"]["score_refs"]
                for j, c in enumerate(bp.categories):
                    if c in out and c in sref:
                        q = _pct(sref[c], dec[j])
                        shipped = {k: out[c].get(k) for k in ("probability", "call", "threshold", "cv_auroc")}
                        out[c] = dict(out[c], probability=round(q, 3),
                                      call="present" if q >= float(out[c]["threshold"]) else "absent", shipped=shipped,
                                      source="bulk_rna (BeatAML-trained), re-referenced to cohort %s" % self.cohort,
                                      cohort_adapter=dict(prov, rereferenced=True))
            zb = {}                              # layer 2: the custom per-mutation models
            for u in self.units:
                if not u["routed"]:
                    continue
                dec = float(u["intercept"])
                for name, coef in u["coef"].items():
                    if name not in zb:
                        zb[name] = self._block_z(name, bp, sample, get)
                    dec += float(np.dot(zb[name], coef))
                self._write(out, u, dec, "custom: " + "+".join(sorted(u["coef"])), prov)
            return out
        xz = np.nan_to_num((xf - self.g_mu) / self.g_sd)
        dz = None
        if any(u["routed"] and u["arm"] == "reweight" for u in self.units):
            dz = (shipped_decisions(bp, sample, self.d["shipped_ref"]) - self.dec_mu) / self.dec_sd
        for u in self.units:
            if u["routed"]:
                feats = dz if u["arm"] == "reweight" else xz
                self._write(out, u, float(np.dot(feats, u["coef"]) + u["intercept"]), u["arm"], prov)
        return out

    def score_units(self, bp, sample):
        """v2 only: {unit name: raw score} for every fitted unit, routed or not (for evaluation)."""
        self._dec_cache = {}
        get = self._linear(bp, sample)
        zb, out = {}, {}
        for u in self.units:
            if "coef" not in u:
                continue
            dec = float(u["intercept"])
            for name, coef in u["coef"].items():
                if name not in zb:
                    zb[name] = self._block_z(name, bp, sample, get)
                dec += float(np.dot(zb[name], coef))
            out[u["name"]] = dec
        return out

    def _write(self, out, u, dec, arm, prov):
        """Place one unit's score into the call dict according to the unit's mode."""
        q = _pct(u["oof_ref"], dec)
        call = "present" if q >= u["threshold"] else "absent"
        row = {"probability": round(q, 3), "call": call, "threshold": round(u["threshold"], 3),
               "cv_auroc": round(u["cv_auroc"], 3), "n_pos_train": u["n_pos"],
               "confidence": "ok" if u["cv_auroc"] >= 0.65 else "abstain: weak (CV AUROC %.2f)" % u["cv_auroc"],
               "source": "cohort adapter %s (%s)" % (self.cohort, arm),
               "cohort_adapter": dict(prov, unit=u["name"], arm=arm, shipped_auroc_on_cohort=u.get("shipped_auroc"))}
        if u["mode"] == "new":
            row["label"] = u.get("label")
            out[u["category"]] = row
            return
        present = [c for c in u["replaces"] if c in out]
        if not present:
            return
        if u["mode"] == "gate":
            # the unit predicts the GENE; it may only veto variant calls of that gene, never assert one
            for c in present:
                note = dict(prov, unit=u["name"], gate=True, gene_probability=round(q, 3), gene_call=call)
                if call == "absent" and out[c].get("call") == "present":
                    shipped = out[c].get("shipped") or {k: out[c].get(k) for k in ("probability", "call", "threshold", "cv_auroc")}
                    out[c] = dict(out[c], call="absent", probability=round(min(out[c].get("probability") or 0, q), 3),
                                  shipped=shipped, source=row["source"], cohort_adapter=dict(note, vetoed=True))
                else:
                    out[c] = dict(out[c], cohort_adapter=note)
            return
        if u["mode"] == "category":
            targets = present
        else:   # gene-level unit over several variant categories: the call goes to the variant the
                # shipped model ranks highest; the adapter predicts the GENE, not which variant
            targets = [max(present, key=lambda c: out[c].get("probability") or 0.0)]
            row["variant_assignment"] = ("gene-level call from the cohort adapter; the variant class is "
                                         "the shipped model's most likely one for this gene")
        for c in present:
            shipped = out[c].get("shipped") or {k: out[c].get(k) for k in ("probability", "call", "threshold", "cv_auroc")}
            if c in targets:
                out[c] = dict(row, shipped=shipped)
            else:
                out[c] = dict(out[c], call="absent", probability=round(min(out[c].get("probability") or 0, q), 3),
                              shipped=shipped, source=row["source"],
                              cohort_adapter=dict(prov, unit=u["name"], superseded_by_gene_call=True))


# ============================================================ dosage block ====
DOSAGE_CLIP, DOSAGE_BURDEN_AT = 2.0, 0.15


def dosage_features(L, d):
    """Expression dosage per chromosome arm for ONE sample.

    L: log2(CP10k+1) over d['genes'] (NaN = not measured). Per gene: log-ratio to the cohort median,
    clipped to +-2, minus the sample's genome-wide median; then the mean over each arm and each whole
    chromosome, plus three burden summaries (sd across arms, largest |arm|, arms beyond +-0.15).
    Monosomy 7 shows as a low chr7 score, trisomy 8 as a high chr8 score, a complex karyotype as burden.
    Gene coordinates were taken from the engine's fastCNV resource when the adapter was built."""
    import warnings
    L = np.asarray(L, float)
    arm, chrom = np.asarray(d["arm"]), np.asarray(d["chrom"])
    R = np.clip(L - np.asarray(d["ref"], float), -DOSAGE_CLIP, DOSAGE_CLIP)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        R = R - np.nanmedian(R)
        A = np.array([np.nanmean(R[arm == a]) for a in d["arms"]])
        C = np.array([np.nanmean(R[chrom == c]) for c in d["chroms"]])
        B = np.array([np.nanstd(A), np.nanmax(np.abs(A)), float((np.abs(A) > DOSAGE_BURDEN_AT).sum())])
    return np.concatenate([A, C, B])


# =============================================================== fitting ====
def fit(bp, X_lin, genes, DEC, labels, units, *, n_splits=5, n_repeats=5, margin=0.03, stability=0.70,
        min_pos=6, new_min_auroc=0.75, n_top=2500, seed=0, n_jobs=1, shipped_calls=None, force_all=False,
        cohort="cohort"):
    """X_lin: patients x genes linear expression (any scale; renormalised to CP10k). DEC: patients x
    categories shipped decision scores (shipped_decisions). labels: {unit_name: array of 0/1/nan}.
    units: list of unit specs. shipped_calls: per-patient shipped call dicts (for the shipped baseline).
    force_all routes every eligible unit to its best arm regardless of the margin (for comparison only).

    Routing is decided on PAIRED inner folds: in every repeat x fold, the cohort model and the shipped
    model are scored on the same held-out training patients. A unit is routed only if the mean fold-wise
    AUROC gain is >= `margin` AND the gain is positive in >= `stability` of the folds AND its F1 is not
    worse. A single pooled comparison was not enough: at ~130 patients it routed units (NRAS) that then
    lost on patients outside the training set."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import StratifiedKFold
    from sklearn.metrics import roc_auc_score

    X = np.vstack([_cp10k(r) for r in np.asarray(X_lin, float)])
    L = np.log1p(X) / LN2
    genes = np.asarray(genes)
    ok = (L.mean(0) > 0.05) & np.array(["." not in g for g in genes])   # unannotated clone lncRNAs are
    cand = np.where(ok)[0]                                                # platform-fragile: excluded
    sel = cand[np.argsort(-L[:, cand].var(0))[:n_top]]
    mu, sd = L[:, sel].mean(0), L[:, sel].std(0)
    sd[sd == 0] = 1.0
    G = (L[:, sel] - mu) / sd
    dmu, dsd = DEC.mean(0), DEC.std(0)
    dsd[dsd == 0] = 1.0
    Dz = (DEC - dmu) / dsd
    mz_train = np.mean(np.abs(G), axis=1)
    pc_train = np.array([np.corrcoef(row, mu)[0, 1] for row in L[:, sel]])

    def model(arm):
        return LogisticRegression(C=0.1 if arm == "reweight" else 1.0, class_weight="balanced", max_iter=3000)

    def f1max(y, s):
        best, thr = 0.0, 1.0
        for t in np.unique(s):
            c = s >= t
            tp = int((c & (y == 1)).sum())
            f = 2 * tp / max(1, int(c.sum()) + int(y.sum()))
            if f > best:
                best, thr = f, float(t)
        return best, thr

    def one(u):
        y_all = np.asarray(labels[u["name"]], float)
        m = np.isfinite(y_all)
        y = y_all[m].astype(int)
        idx = np.where(m)[0]
        spec = {k: u[k] for k in u}
        spec.update({"n_pos": int(y.sum()), "n": int(len(y)), "routed": False})
        if y.sum() < min_pos or (len(y) - y.sum()) < min_pos:
            spec["reason"] = "too few labelled positives/negatives in training (%d/%d)" % (y.sum(), len(y) - y.sum())
            return spec
        # shipped baseline on the same patients
        sh_auc = sh_f1 = sp = None
        if u["mode"] != "new" and shipped_calls is not None:
            cats = u["replaces"]
            sp = np.array([max((shipped_calls[i].get(c) or {}).get("probability") or 0.0 for c in cats) for i in idx])
            sh_auc = float(roc_auc_score(y, sp))
            sh_f1 = f1max(y, sp)[0]
        res = {}
        for arm in ("reweight", "scratch"):
            Xa = Dz[idx] if arm == "reweight" else G[idx]
            oof = np.zeros((n_repeats, len(y)))
            aucs, fold_auc, fold_diff = [], [], []
            for r in range(n_repeats):
                for tr, te in StratifiedKFold(n_splits, shuffle=True, random_state=seed + r).split(Xa, y):
                    oof[r, te] = model(arm).fit(Xa[tr], y[tr]).decision_function(Xa[te])
                    if 0 < y[te].sum() < len(te):
                        a_c = roc_auc_score(y[te], oof[r, te])
                        fold_auc.append(a_c)
                        if sp is not None:
                            fold_diff.append(a_c - roc_auc_score(y[te], sp[te]))
                aucs.append(roc_auc_score(y, oof[r]))
            # percentile of each patient's OOF score within its repeat, then averaged across repeats
            q = np.mean([np.searchsorted(np.sort(o), o, side="right") / len(o) for o in oof], axis=0)
            f1, thr = f1max(y, q)
            res[arm] = {"auc": float(np.mean(aucs)), "f1": f1, "thr": thr, "oof": oof.mean(0),
                        "fold_auc": np.asarray(fold_auc), "fold_diff": np.asarray(fold_diff)}
        arm = max(res, key=lambda a: res[a]["auc"])
        best = res[arm]
        spec.update({"arm": arm, "cv_auroc": round(best["auc"], 4), "cv_f1": round(best["f1"], 4),
                     "shipped_auroc": None if sh_auc is None else round(sh_auc, 4),
                     "shipped_f1_opt": None if sh_f1 is None else round(sh_f1, 4),
                     "arms": {a: {"cv_auroc": round(v["auc"], 4), "cv_f1": round(v["f1"], 4)} for a, v in res.items()}})
        if u["mode"] == "new":
            stab = float(np.mean(best["fold_auc"] >= 0.6)) if best["fold_auc"].size else 0.0
            go = best["auc"] >= new_min_auroc and stab >= stability
            spec["fold_stability"] = round(stab, 3)
            why = ("%s: cohort CV AUROC %.3f (floor %.2f), >=0.6 in %.0f%% of folds (need %.0f%%)"
                   % ("added" if go else "not added", best["auc"], new_min_auroc, 100 * stab, 100 * stability))
        else:
            d = best["fold_diff"]
            gain = float(d.mean()) if d.size else 0.0
            stab = float(np.mean(d > 0)) if d.size else 0.0
            spec.update({"fold_gain": round(gain, 4), "fold_stability": round(stab, 3)})
            go = force_all or (gain >= margin and stab >= stability and best["f1"] >= sh_f1)
            why = ("%s: fold-paired AUROC gain %+.3f (need +%.2f), better in %.0f%% of %d folds (need %.0f%%), "
                   "F1 %.3f vs shipped %.3f"
                   % ("routed" if go else "kept shipped", gain, margin, 100 * stab, d.size, 100 * stability,
                      best["f1"], sh_f1))
        spec["reason"] = why
        if not go:
            return spec
        Xa = Dz[idx] if arm == "reweight" else G[idx]
        mdl = model(arm).fit(Xa, y)
        spec.update({"routed": True, "coef": [float(c) for c in mdl.coef_[0]],
                     "intercept": float(mdl.intercept_[0]), "threshold": float(best["thr"]),
                     "oof_ref": sorted(float(v) for v in best["oof"])})
        return spec

    if n_jobs and n_jobs != 1:
        from joblib import Parallel, delayed
        specs = Parallel(n_jobs=n_jobs)(delayed(one)(u) for u in units)
    else:
        specs = [one(u) for u in units]
    d = {"version": VERSION, "cohort": cohort, "shipped_ref": "sc",
         "selection": {"n_splits": n_splits, "n_repeats": n_repeats, "margin_auroc": margin,
                       "stability": stability,
                       "rule": "route iff the fold-paired AUROC gain over shipped is >= margin on average AND "
                               "positive in >= stability of the folds AND the F1 (both at their F1-optimal cut "
                               "on training data) is not worse; add a new unit iff CV AUROC >= floor and "
                               "fold AUROC >= 0.6 in >= stability of the folds",
                       "min_pos": min_pos, "new_min_auroc": new_min_auroc, "force_all": bool(force_all)},
         "features": {"genes": [str(g) for g in genes[sel]], "mu": mu.tolist(), "sd": sd.tolist(),
                      "dec_mu": dmu.tolist(), "dec_sd": dsd.tolist(), "n_train": int(len(X))},
         "guards": {"min_gene_coverage": 0.90,
                    "max_mean_abs_z": float(max(2.0, 1.5 * np.quantile(mz_train, 0.99))),
                    # 5th percentile of the training samples, never looser than 0.85: refuses a few
                    # low-quality in-cohort samples (they keep the shipped calls) rather than risk
                    # rewriting calls for another platform
                    "min_profile_corr": float(min(0.85, np.quantile(pc_train, 0.05)))},
         "units": specs}
    return CohortAdapter(d)
