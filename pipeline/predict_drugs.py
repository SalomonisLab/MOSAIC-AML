#!/usr/bin/env python3
"""COMPASS-AML patient entry point: expression in, tiered drug prioritisation out.

Runs the whole layer for one sample and writes a report:

  Model A   response prediction for every modelled inhibitor, calibrated and percentiled
  Model B   the same model re-applied per cell state (single-cell input only) -> coverage, escape
  Model C   mechanistic target/pathway evidence, computed independently of A
  utility   S_ij, with every penalty term itemised
  agents    the eight expert agents, all non-voting with respect to the anchored subtype call

Accepts either input the platform already produces:
  --pseudobulk  a (cell-state x gene) raw-count TSV, or an atlas sample name to pull from the h5ad
  --bulk        a single expression vector (TSV: gene <tab> value)

  python predict_drugs.py --atlas-sample "3v2::GSM3901485" -o runs/<id>/drug_report.json
"""
import os, sys, json, time, pickle, argparse, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd

from amlmm.drug import (data as D, model as M, statemodel as SM, mechanism as MC,
                        targets as TG, utility as U, agents as AG)

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
MODEL = os.path.join(HERE, "drug_response_model.pkl")
BUNDLE = os.path.join(os.path.dirname(ROOT), "aml-bakeoff", "bundle_data.npz")
H5 = os.path.join(ROOT, "data", "RNA", "pseudobulk_counts_hashed.h5ad")
STATE_COL = "Hs-BM-titrated-reference-centroid"


# ------------------------------------------------------------------ input ----
def atlas_sample(name):
    """Pull one atlas sample's (cell-state x gene) raw counts straight from the pseudobulk h5ad,
    reading only that sample's rows rather than materialising the 1 GB matrix."""
    import h5py
    from amlmm.drug.h5rows import obs_column, var_index, read_rows
    with h5py.File(H5, "r") as f:
        ds, sm = obs_column(f, "Dataset"), obs_column(f, "Sample")
        key = np.array(["%s::%s" % (a, b) for a, b in zip(ds, sm)])
        m = np.where(key == name)[0]
        if not len(m):
            raise SystemExit("sample %r not in the atlas (e.g. %s)" % (name, sorted(set(key))[:3]))
        st = obs_column(f, STATE_COL)[m]
        nc = np.asarray(obs_column(f, "n_cells"), dtype=float)[m]
        genes = var_index(f)
    _, X = read_rows(H5, m)
    counts = pd.DataFrame(X, index=st, columns=genes).groupby(level=0).sum()
    return counts, pd.Series(nc, index=st).groupby(level=0).sum()


# ------------------------------------------------------------------- run ----
def run(state_counts=None, n_cells=None, bulk=None, mutations=None, clinical=None,
        model_path=MODEL, top=10, prob_cut=0.5):
    with open(model_path, "rb") as f:
        mod = pickle.load(f)
    fs = mod.fs
    drugs = list(mod.drug_models)
    d = np.load(BUNDLE, allow_pickle=True)

    # ---- the whole-sample expression vector, however it arrived ----
    if state_counts is not None:
        lin = SM.cp10k(state_counts.values.sum(0, keepdims=True))
        vec = pd.DataFrame(lin, index=["bulk"], columns=state_counts.columns)
    elif bulk is not None:
        vec = pd.DataFrame([bulk.values], index=["bulk"], columns=bulk.index)
    else:
        raise SystemExit("need --pseudobulk / --atlas-sample or --bulk")
    # Which cohort is this sample? Everything downstream -- the expression z-reference, the score
    # percentile reference, the mechanistic percentiles, the OOD distance -- must be matched to it.
    is_sc = state_counts is not None
    expr_ref = "sc" if (is_sc and "sc" in fs.ref) else "beataml"
    score_ref = "sc_sample" if (is_sc and "sc_sample" in (getattr(mod, "score_refs", {}) or {})) else "beataml"

    sr = SM.StateResponse(mod)
    Zb, n_shared = sr._z(vec, ref=expr_ref)

    # ---- Model A ----
    preds = mod.predict_patient(Zb, mut=None, meta=None, drugs=drugs, ref=score_ref)

    # ---- Model B ----
    state_result = None
    if is_sc:
        state_result = sr.predict(state_counts, n_cells, drugs=drugs, prob_cut=prob_cut)

    # ---- Model C ----
    mech_ref = d["sc_X"].astype(float) if is_sc else d["ba_X"].astype(float)
    mech_model = MC.MechanismModel(fs.genes, mech_ref, sym2ens=mod.sym2ens)
    x_lin = np.zeros(len(fs.genes))
    gidx = {g: i for i, g in enumerate(fs.genes)}
    for g, v in zip(vec.columns, vec.values[0]):
        g = str(g); k = g if g in gidx else mod.sym2ens.get(g)
        if k in gidx:
            x_lin[gidx[k]] = v
    state_scores = None
    if state_result:
        tot = sum(s["fraction"] for s in state_result["states"]) or 1.0
        state_scores = {}
        for s in state_result["states"]:
            state_scores[s["group"]] = state_scores.get(s["group"], 0.0) + s["fraction"] / tot
    mech = {dr: mech_model.evaluate(dr, x_lin, mutations, state_scores) for dr in drugs}

    # ---- OOD + differentiation-state percentile, both against the MATCHED cohort ----
    nb = mod.neighbours(Zb, k=20, cohort=("sc" if is_sc else "beataml"))
    nnb = getattr(mod, "nn", None) or {}
    ood_ref = nnb.get("ood_ref_sc") if (is_sc and "ood_ref_sc" in nnb) else nnb.get("ood_ref")
    ood_dist = None if nb is None else nb["self_distance"]
    ood_q = None if (ood_ref is None or ood_dist is None) else float((ood_ref <= ood_dist).mean())
    Sb, snames = fs._state_block(Zb)
    axis = float(Sb[0, snames.index("axis_prim_minus_mature")])
    Zref = fs.z(mech_ref, expr_ref)
    Sref, _ = fs._state_block(Zref)
    axis_q = float((Sref[:, snames.index("axis_prim_minus_mature")] <= axis).mean())

    # ---- utility + agents ----
    oofm = getattr(mod, "oof_metrics", {}) or {}
    tiers = getattr(mod, "drug_tier", {}) or {}
    scored, abstained, per_ev = {}, [], {}
    resp_agent = AG.drug_response_agent(mod, preds, nb)
    resp_by_drug = {r["inhibitor"]: r for r in resp_agent["evidence"]["all"]}
    for dr in drugs:
        p = preds[dr]
        sm_ = (state_result or {}).get("per_drug", {}).get(dr)
        s = U.score(dr, p.get("prob_sensitive"), (oofm.get(dr) or {}).get("auroc"), p.get("n_train"),
                    mech_evidence=mech.get(dr), state_metrics=sm_, clinical=clinical,
                    ood_distance=ood_dist, ood_reference=ood_ref)
        ch = AG.skeptic_agent(mod, dr, p, mech.get(dr), sm_, state_result, ood_q,
                              tiers.get(dr), axis_q, curves=None)
        s["challenges"] = ch
        if s["components"]["uncertainty"] >= 0.75:
            abstained.append({"inhibitor": dr, "reason": "uncertainty %.2f above the abstention "
                                                         "threshold" % s["components"]["uncertainty"]})
            continue
        scored[dr] = s
        per_ev[dr] = {"response": resp_by_drug.get(dr), "state": sm_, "mechanism": mech.get(dr),
                      "challenges": ch}

    ranked = U.rank_by_tier(scored)
    agents = [resp_agent,
              AG.cell_state_agent(state_result),
              AG.mechanism_agent(mech),
              AG.clinical_evidence_agent(drugs, clinical),
              AG.combination_agent(state_result, scored),
              AG.reporting_agent(ranked, per_ev, abstained,
                                 {"genes_shared_with_model": int(n_shared),
                                  "ood_percentile": None if ood_q is None else round(ood_q, 3),
                                  "differentiation_axis_percentile": round(axis_q, 3),
                                  "nearest_beataml_specimens": None if nb is None else nb["specimens"][:10]},
                                 top_per_tier=top)]
    return {"generated": time.strftime("%Y-%m-%d %H:%M"),
            "model_version": M.DrugResponseModel.VERSION,
            "n_drugs_modelled": len(drugs), "n_drugs_reported": len(scored),
            "n_abstained": len(abstained),
            "patient": {"genes_shared_with_model": int(n_shared),
                        "ood_distance": ood_dist, "ood_percentile": ood_q,
                        "differentiation_axis_percentile": round(axis_q, 3),
                        "single_cell_states": None if state_result is None else len(state_result["states"])},
            "ranked": ranked, "abstained": abstained, "agents": agents,
            "per_drug": per_ev}


def markdown(rep, title="COMPASS-AML drug prioritisation"):
    L = ["# %s" % title, "",
         "*Predicted **ex-vivo** sensitivity from the BeatAML2 functional screen. This is a "
         "prioritisation for trial matching or laboratory validation, **not** a treatment "
         "recommendation.*", ""]
    p = rep["patient"]
    L += ["- inhibitors modelled: **%d**, reported %d, abstained %d"
          % (rep["n_drugs_modelled"], rep["n_drugs_reported"], rep["n_abstained"]),
          "- genes shared with the model: %d" % p["genes_shared_with_model"],
          "- distance from the BeatAML training distribution: %s percentile"
          % ("%.0fth" % (100 * p["ood_percentile"]) if p["ood_percentile"] is not None else "n/a"),
          "- differentiation axis (primitive - mature): %.0fth percentile of BeatAML"
          % (100 * p["differentiation_axis_percentile"]),
          "- cell states scored: %s" % (p["single_cell_states"] or "bulk input, none"), ""]
    for tier, blk in rep["ranked"].items():
        if not blk["ranked"]:
            continue
        L += ["## %s  (%d considered)" % (blk["label"], blk["n"]), "",
              "| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |",
              "|---|---|---|---|---|---|---|---|---|"]
        for i, r in enumerate(blk["ranked"][:10], 1):
            c = r["components"]
            L.append("| %d | %s | %.3f | %s | %s | %s | %s | %s | %d |"
                     % (i, r["inhibitor"], r["utility"],
                        _f(c["sensitivity"]), _f((rep["per_drug"].get(r["inhibitor"], {}).get("response") or {}).get("model_oof_auroc")),
                        _f(c["coverage_blast"]), _f(c["coverage_LSC_like"]), _f(c["mechanistic"]),
                        len(r.get("challenges") or [])))
        L.append("")
    comb = next((a for a in rep["agents"] if a["name"] == "combination"), None)
    if comb and comb["evidence"].get("pairs"):
        L += ["## Combination hypotheses (complementary cell-state coverage)", "",
              "*BeatAML2 measured single agents only; no synergy is claimed.*", "",
              "| pair | pathways | coverage A | coverage B | union | gain |", "|---|---|---|---|---|---|"]
        for c in comb["evidence"]["pairs"]:
            L.append("| %s + %s | %s / %s | %.2f | %.2f | %.2f | +%.2f |"
                     % (c["pair"][0], c["pair"][1], c["pathways"][0], c["pathways"][1],
                        c["coverage_a"], c["coverage_b"], c["coverage_union"], c["gain"]))
        L.append("")
    return "\n".join(L)


def _f(v):
    return "—" if v is None else ("%.3f" % v)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--atlas-sample")
    ap.add_argument("--pseudobulk", help="(cell-state x gene) raw-count TSV")
    ap.add_argument("--bulk", help="two-column TSV: gene <tab> value")
    ap.add_argument("--mutations", help="JSON {category: probability}")
    ap.add_argument("--clinical", help="JSON {age:.., prior_lines:..}")
    ap.add_argument("-o", "--out", default=os.path.join(ROOT, "deliverables", "drug_report_example.json"))
    a = ap.parse_args()

    sc = nc = bulk = None
    if a.atlas_sample:
        sc, nc = atlas_sample(a.atlas_sample)
    elif a.pseudobulk:
        sc = pd.read_csv(a.pseudobulk, sep="\t", index_col=0)
        nc = pd.Series(1000.0, index=sc.index)
    elif a.bulk:
        b = pd.read_csv(a.bulk, sep="\t", index_col=0)
        bulk = b.iloc[:, 0]
    mut = json.load(open(a.mutations)) if a.mutations else None
    clin = json.load(open(a.clinical)) if a.clinical else None

    rep = run(sc, nc, bulk, mut, clin)
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    json.dump(rep, open(a.out, "w"), indent=1, default=str)
    md = a.out.replace(".json", ".md")
    open(md, "w", encoding="utf-8").write(markdown(rep))
    print(markdown(rep))
    print("\nwrote %s + %s" % (a.out, md))


if __name__ == "__main__":
    main()
