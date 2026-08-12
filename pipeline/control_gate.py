#!/usr/bin/env python3
"""Deployed HEALTHY-vs-DISEASED control gate.

A new scRNA upload yields cell-state COMPOSITION (cellHarmony) directly — the imputed modalities
(Lipid/Metabolite/GRN/ADT) are not wired for uploads — and composition separates healthy from diseased
at AUC ~0.91 within-dataset (validated, biology not batch). This gate freezes a composition->control/diseased
classifier and is applied to a new sample BEFORE mutation calling: a sample called 'control' is reported
as "control / no mutation".

Operating point is CONSERVATIVE by design: the threshold keeps disease SENSITIVITY >= 0.95 (never skip a
real patient), accepting that some true controls are called 'diseased' — those simply come back clean from
the mutation panel. So a 'control' call is high-confidence healthy.

  train:  python control_gate.py --train         (on an LSF compute node)
  use:    from control_gate import load_gate, score_gate
-> pipeline/control_gate.pkl
"""
import os, sys, json, pickle, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd

GATE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "control_gate.pkl")
MODALITY = "Composition"
TARGET_DISEASE_SENS = 0.99                                # never miss a real patient; mis-gated controls just proceed and come back clean


def _labels(ctx):
    s = ctx.tables["samples"]; dc = s.get("disease_category").astype(str)
    y = pd.Series(index=s.index, dtype=float)
    y[dc.eq("Control")] = 0.0
    y[dc.isin({"AML", "MDS", "T-ALL"})] = 1.0
    return y.dropna()


def _composition(ctx):
    from amlmm import discovery as D
    cache = ctx.path("_sl_Composition.pkl")
    B = pd.read_pickle(cache) if os.path.exists(cache) else D._sample_level_matrix(ctx, "composition", set(ctx.tables["samples"].index))
    return B.fillna(0.0)


def _thresh_for_sensitivity(y, p, target):
    """Largest p_diseased cutoff t (predict diseased if p>=t) that still keeps disease sensitivity >= target."""
    best = 0.0
    for t in np.unique(p):
        sens = ((p >= t) & (y == 1)).sum() / max(1, (y == 1).sum())
        if sens >= target:
            best = float(t)
    return best


def train(ctx):
    from sklearn.preprocessing import StandardScaler
    from sklearn.model_selection import GroupKFold
    from sklearn.metrics import roc_auc_score
    from sklearn.linear_model import LogisticRegression
    y = _labels(ctx); B = _composition(ctx); dg = ctx.tables["samples"]["donor_group"].astype(str)
    ids = [s for s in B.index if s in y.index]
    X = B.loc[ids].values; cols = list(B.columns); yy = y.loc[ids].values.astype(int); g = dg.loc[ids].values

    oof = np.full(len(ids), np.nan)                                   # donor-grouped CV -> honest AUC + threshold
    for tri, vai in GroupKFold(5).split(X, yy, g):
        if len(set(yy[tri])) < 2:
            continue
        keep = X[tri].std(0) > 0
        sc = StandardScaler().fit(X[tri][:, keep])
        m = LogisticRegression(C=0.1, class_weight="balanced", max_iter=3000).fit(sc.transform(X[tri][:, keep]), yy[tri])
        oof[vai] = m.predict_proba(sc.transform(X[vai][:, keep]))[:, 1]
    ok = ~np.isnan(oof); p, yt = oof[ok], yy[ok]
    auc = roc_auc_score(yt, p)
    thr = _thresh_for_sensitivity(yt, p, TARGET_DISEASE_SENS)
    sens = ((p >= thr) & (yt == 1)).sum() / max(1, (yt == 1).sum())
    spec = ((p < thr) & (yt == 0)).sum() / max(1, (yt == 0).sum())

    keep = X.std(0) > 0                                               # final model on the FULL cohort
    sc = StandardScaler().fit(X[:, keep])
    model = LogisticRegression(C=0.1, class_weight="balanced", max_iter=3000).fit(sc.transform(X[:, keep]), yy)
    # Keep the out-of-fold scores. Without them any question about the operating point ("what would a
    # balanced threshold give?") needed a full retrain with atlas access, so the 0.50 control
    # specificity went unexamined for months while the gate's positive call was labelled 'diseased'.
    j = (np.array([((p >= t) & (yt == 1)).sum() / max(1, (yt == 1).sum())
                   + ((p < t) & (yt == 0)).sum() / max(1, (yt == 0).sum()) - 1 for t in np.unique(p)]))
    t_youden = float(np.unique(p)[int(np.argmax(j))]) if len(j) else float(thr)
    alt = {}
    for nm, t in (("youden", t_youden), ("balanced_0.5", 0.5)):
        alt[nm] = {"threshold": round(float(t), 4),
                   "sens_disease": round(float(((p >= t) & (yt == 1)).sum() / max(1, (yt == 1).sum())), 4),
                   "spec_control": round(float(((p < t) & (yt == 0)).sum() / max(1, (yt == 0).sum())), 4)}
    gate = {"modality": MODALITY, "cols": cols, "keep": keep, "scaler": sc, "model": model,
            "threshold": float(thr), "cv_auc": float(auc), "cv_sens_disease": float(sens),
            "cv_spec_control": float(spec), "n": len(ids), "n_control": int((yy == 0).sum()),
            "oof_scores": p.tolist(), "oof_labels": yt.tolist(),
            "alternative_operating_points": alt}
    print("  alternative operating points: %s" % json.dumps(alt))
    pickle.dump(gate, open(GATE_PATH, "wb"))
    print("CONTROL GATE trained: cv_auc=%.3f thr=%.3f sens(disease)=%.3f spec(control)=%.3f n=%d nControl=%d"
          % (auc, thr, sens, spec, len(ids), int((yy == 0).sum())))
    # sanity: how the frozen gate calls the cohort controls vs diseased (in-sample, illustrative)
    pall = model.predict_proba(sc.transform(X[:, keep]))[:, 1]
    call = (pall >= thr).astype(int)
    print("  cohort controls called 'control': %d/%d ; diseased called 'diseased': %d/%d -> %s"
          % (int(((call == 0) & (yy == 0)).sum()), int((yy == 0).sum()),
             int(((call == 1) & (yy == 1)).sum()), int((yy == 1).sum()), GATE_PATH))
    return gate


def load_gate(path=GATE_PATH):
    return pickle.load(open(path, "rb")) if os.path.exists(path) else None


def score_gate(gate, comp_series):
    """comp_series: cell-state -> fraction for ONE sample.

    The positive call is `not_excluded`, NOT `diseased`, and the difference is the whole point. The
    threshold is chosen to keep disease sensitivity at 0.99, and the price of that is a measured
    **control specificity of 0.50** on 24 training controls: the gate calls half of genuinely healthy
    specimens positive. Two normal sorted populations in GSE281087 (CD34-1, GMP-1) were called
    'diseased' by the old label, which is not a malfunction -- it is what a 0.50-specificity gate does,
    and reporting it as 'diseased' overstated it every time.

    So: a **control** call is strong evidence of health and is worth acting on (that direction is
    99%-sensitivity-protected). A **not_excluded** call is close to no evidence at all, and every field
    below exists so a reader cannot mistake it for a positive finding.
    """
    cols = gate["cols"]
    vec = np.array([float(comp_series.get(c.split("::", 1)[-1], comp_series.get(c, 0.0))) for c in cols], float)
    s = vec.sum()
    if s > 0:
        vec = vec / s
    Xk = vec[gate["keep"]].reshape(1, -1)
    p = float(gate["model"].predict_proba(gate["scaler"].transform(Xk))[:, 1][0])
    pos = p >= gate["threshold"]
    spec = gate.get("cv_spec_control")
    return {"call": "not_excluded" if pos else "control",
            "p_diseased": round(p, 4), "threshold": round(gate["threshold"], 4),
            "cv_auc": gate["cv_auc"], "modality": gate["modality"],
            "control_specificity": spec, "n_control_train": gate.get("n_control"),
            "interpretation": (
                ("disease NOT excluded. This is not a positive finding: at the deployed operating "
                 "point the gate calls %.0f%% of genuinely healthy specimens positive, so this call "
                 "carries little information on its own." % (100 * (1 - (spec if spec is not None else 0.5))))
                if pos else
                ("called healthy. The threshold is set for %.0f%% disease sensitivity, so a control "
                 "call is high-confidence." % (100 * (gate.get("cv_sens_disease") or 0.99)))),
            "trained_on": ("%s samples, of which %s controls -- the control side is the thin one"
                           % (gate.get("n"), gate.get("n_control")))}


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--train", action="store_true"); a = ap.parse_args()
    from amlmm.context import build_context, Config
    ctx = build_context(Config(run_id="single_modality"))
    if a.train:
        train(ctx)
    else:
        print("use --train")
