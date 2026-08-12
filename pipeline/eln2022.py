#!/usr/bin/env python3
"""ELN 2022 risk classification, inferred for BeatAML2.

BeatAML ships an `ELN2017` column and nothing later, so every benchmark in this platform has been
scored against the 2017 standard while claiming relevance to current practice. The two differ in ways
that matter for exactly the patients we get wrong:

  FLT3-ITD allelic ratio   USED in 2017 (AR >= 0.5 changed the category), DROPPED in 2022. Any
                           FLT3-ITD with mutated NPM1 is intermediate, ratio irrelevant.
  CEBPA                    2017 required BIALLELIC. 2022 requires an in-frame mutation in the bZIP
                           region, monoallelic or biallelic -- a different set of patients.
  MDS-related genes        NEW in 2022: ASXL1, BCOR, EZH2, RUNX1, SF3B1, SRSF2, STAG2, U2AF1, ZRSR2
                           are adverse (2017 had only ASXL1 and RUNX1).
  TP53                     NEW in 2022, and explicitly AT A VARIANT ALLELE FRACTION >= 10%,
                           irrespective of allelic status.
  t(9;11)                  intermediate, and takes precedence over concurrent adverse gene mutations.

Source: Dohner H, Wei AH, Appelbaum FR, et al. Diagnosis and management of AML in adults: 2022 ELN
recommendations. Blood 140(12):1345-1377, Table 6.

THE VAF THRESHOLD IS A REAL DEGREE OF FREEDOM. The guideline names 10% for TP53 and is silent for
every other gene. 40% is the conventional proxy for a clonal/biallelic event. In BeatAML, TP53 is
mutated in 80 specimens at >= 10% and 70 at >= 40%, so the choice moves one in eight TP53 patients --
which is why `--vaf` is a parameter here and every downstream benchmark is run at both.

  python eln2022.py [--vaf 0.10] [--out ...]  ->  labels/eln2022_beataml_vaf10.tsv
"""
import os, sys, re, json, argparse, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
BA = os.path.join(ROOT, "data", "external", "beataml")

MR_GENES = ["ASXL1", "BCOR", "EZH2", "RUNX1", "SF3B1", "SRSF2", "STAG2", "U2AF1", "ZRSR2"]
# CEBPA is 358 aa; the basic leucine zipper spans roughly 278-358 (basic region ~278-296, zipper to 358)
CEBPA_BZIP = (278, 358)
INFRAME = {"inframe_insertion", "inframe_deletion", "missense_variant"}


# ---------------------------------------------------------------- karyotype parsing
def _clones(k):
    """Split a karyotype string into clone strings, dropping the cell counts in [ ]."""
    k = re.sub(r"\[[^\]]*\]", "", str(k))
    return [c.strip() for c in re.split(r"/", k) if c.strip()]


def karyo_flags(k):
    """Structured cytogenetic flags from a free-text ISCN karyotype."""
    f = dict.fromkeys(("normal", "t_8_21", "inv16", "t_9_11", "t_6_9", "t_9_22", "t_8_16",
                       "inv3", "mecom", "kmt2a_r", "minus5_del5q", "minus7", "abn17p",
                       "complex", "monosomal", "unknown"), False)
    if not isinstance(k, str) or not k.strip() or k.strip().lower() in ("na", "nan", "unknown"):
        f["unknown"] = True
        return f
    s = k.replace(" ", "")
    f["normal"] = bool(re.search(r"^4[56],X[XY]$", _clones(s)[0] if _clones(s) else ""))
    f["t_8_21"] = bool(re.search(r"t\(8;21\)", s))
    f["inv16"] = bool(re.search(r"inv\(16\)|t\(16;16\)", s))
    f["t_9_11"] = bool(re.search(r"t\(9;11\)", s))
    f["t_6_9"] = bool(re.search(r"t\(6;9\)", s))
    f["t_9_22"] = bool(re.search(r"t\(9;22\)", s))
    f["t_8_16"] = bool(re.search(r"t\(8;16\)", s))
    f["inv3"] = bool(re.search(r"inv\(3\)|t\(3;3\)", s))
    f["mecom"] = bool(re.search(r"3q26", s))
    # KMT2A-rearranged = any 11q23 translocation. t(9;11) is handled separately (intermediate).
    f["kmt2a_r"] = bool(re.search(r"t\([^)]*;11\)\([^)]*q23", s) or re.search(r"11q23", s))
    f["minus5_del5q"] = bool(re.search(r"-5\b|del\(5\)\(q|del\(5q", s))
    f["minus7"] = bool(re.search(r"-7\b(?!p)|del\(7\)\(q|monosomy7", s))
    f["abn17p"] = bool(re.search(r"-17\b|del\(17\)\(p|i\(17\)\(q|add\(17\)\(p|17p", s))

    # complex: >=3 unrelated abnormalities in the largest clone, excluding pure hyperdiploidy
    best = 0
    for c in _clones(s):
        body = re.sub(r"^\d+,X[XY]?,?", "", c)
        parts = [p for p in re.split(r",(?![^()]*\))", body) if p]
        n = len(parts)
        structural = any(re.search(r"t\(|inv\(|del\(|add\(|der\(|i\(|dup\(|ins\(", p) for p in parts)
        trisomies = sum(1 for p in parts if re.match(r"^\+\d+$", p))
        if trisomies >= 3 and not structural:
            continue                                  # hyperdiploid without structural change: excluded
        best = max(best, n)
    f["complex"] = best >= 3

    # monosomal: >=2 autosomal monosomies, or 1 autosomal monosomy + >=1 structural abnormality
    for c in _clones(s):
        mono = re.findall(r"-(\d+)\b", c)
        mono = [x for x in mono if x not in ("X", "Y")]
        structural = bool(re.search(r"t\(|inv\(|del\(|add\(|der\(|i\(|dup\(|ins\(", c))
        if len(mono) >= 2 or (len(mono) == 1 and structural):
            f["monosomal"] = True
    return f


# ---------------------------------------------------------------- per-patient classification
def classify(row, muts, vaf):
    """ELN 2022 category for one specimen. Returns (risk, reasons)."""
    why = []
    k = row["_karyo"]
    fus = str(row.get("consensusAMLFusions") or "")
    npm1 = str(row.get("NPM1", "")).lower().startswith("pos")
    itd = str(row.get("FLT3-ITD", "")).lower().startswith("pos")

    def fus_has(*pats):
        return any(re.search(p, fus, re.I) for p in pats)

    # --- class-defining favorable cytogenetics (KIT/FLT3 co-mutation does not alter these) ---
    cbf = k["t_8_21"] or k["inv16"] or fus_has(r"RUNX1[-:]+RUNX1T1", r"CBFB[-:]+MYH11")
    # --- adverse cytogenetics -------------------------------------------------------------
    adv_cyto = []
    if k["t_6_9"] or fus_has(r"DEK[-:]+NUP214"): adv_cyto.append("t(6;9)/DEK::NUP214")
    if k["t_9_22"] or fus_has(r"BCR[-:]+ABL1"): adv_cyto.append("t(9;22)/BCR::ABL1")
    if k["t_8_16"] or fus_has(r"KAT6A[-:]+CREBBP"): adv_cyto.append("t(8;16)/KAT6A::CREBBP")
    if k["inv3"] or k["mecom"] or fus_has(r"MECOM", r"EVI1"): adv_cyto.append("inv(3)/MECOM")
    if k["kmt2a_r"] and not k["t_9_11"]: adv_cyto.append("KMT2A-rearranged")
    if k["minus5_del5q"]: adv_cyto.append("-5/del(5q)")
    if k["minus7"]: adv_cyto.append("-7")
    if k["abn17p"]: adv_cyto.append("-17/abn(17p)")
    if k["complex"]: adv_cyto.append("complex karyotype")
    if k["monosomal"]: adv_cyto.append("monosomal karyotype")

    # --- molecular ------------------------------------------------------------------------
    g = muts.get(row["_spec"], {})
    tp53 = any(v >= vaf for v in g.get("TP53", []))
    mr_hit = sorted([x for x in MR_GENES if any(v >= vaf for v in g.get(x, []))])
    cebpa_bzip = g.get("_cebpa_bzip", False)

    # ---------------- rule order (Table 6 footnotes) --------------------------------------
    # t(9;11) takes precedence over rare concurrent adverse-risk gene mutations
    if k["t_9_11"] or fus_has(r"MLLT3[-:]+KMT2A"):
        why.append("t(9;11)/MLLT3::KMT2A takes precedence over concurrent adverse mutations")
        return "Intermediate", why

    if cbf:
        why.append("core-binding-factor AML (KIT/FLT3 co-mutation does not alter the category)")
        if mr_hit:
            why.append("MDS-related mutations (%s) NOT counted adverse alongside a favorable subtype"
                       % ",".join(mr_hit))
        return "Favorable", why

    if npm1 and adv_cyto:
        why.append("mutated NPM1 with adverse-risk cytogenetics (%s) -> adverse" % "; ".join(adv_cyto))
        return "Adverse", why
    if npm1 and not itd:
        why.append("mutated NPM1 without FLT3-ITD")
        if mr_hit:
            why.append("MDS-related mutations (%s) NOT counted adverse alongside a favorable subtype"
                       % ",".join(mr_hit))
        return "Favorable", why
    if cebpa_bzip and not adv_cyto:
        why.append("in-frame bZIP CEBPA mutation (monoallelic or biallelic, per ELN 2022)")
        return "Favorable", why

    if tp53:
        why.append("TP53 mutated at VAF >= %.0f%%" % (100 * vaf))
        return "Adverse", why
    if adv_cyto:
        why.append("adverse cytogenetics: " + "; ".join(adv_cyto))
        return "Adverse", why
    if mr_hit:
        why.append("MDS-related gene mutation: " + ",".join(mr_hit))
        return "Adverse", why

    if npm1 and itd:
        why.append("mutated NPM1 with FLT3-ITD (allelic ratio not used in ELN 2022)")
        return "Intermediate", why
    if itd:
        why.append("wild-type NPM1 with FLT3-ITD, no adverse lesion")
        return "Intermediate", why
    why.append("no favorable or adverse defining lesion")
    return "Intermediate", why


def load_mutations(vaf_floor=0.0):
    m = pd.read_csv(os.path.join(BA, "mutations.txt"), sep="\t", low_memory=False)
    m["t_vaf"] = pd.to_numeric(m["t_vaf"], errors="coerce")
    m = m[m["t_vaf"].notna()]
    out = {}
    for spec, sub in m.groupby("dbgap_sample_id"):
        d = {}
        for sym, s2 in sub.groupby("symbol"):
            d[str(sym)] = list(s2["t_vaf"].values)
        # CEBPA bZIP, in-frame only, position inside the bZIP window
        ce = sub[sub["symbol"] == "CEBPA"]
        bz = False
        for _, r in ce.iterrows():
            if str(r["variant_classification"]) not in INFRAME:
                continue
            pp = str(r.get("protein_position") or "")
            nums = [int(x) for x in re.findall(r"\d+", pp.split("/")[0])]
            if nums and CEBPA_BZIP[0] <= max(nums) <= CEBPA_BZIP[1] and r["t_vaf"] >= vaf_floor:
                bz = True
        d["_cebpa_bzip"] = bz
        out[str(spec)] = d
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vaf", type=float, default=0.10)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    cl = pd.read_excel(os.path.join(BA, "beataml_wv1to4_clinical.xlsx"))
    cl["_spec"] = cl["dbgap_rnaseq_sample"].astype(str)
    # The clinical file has 942 rows but only 699 distinct RNA-seq specimens (repeats and a shared
    # 'nan' key). Reporting over the raw rows double-counts patients and shifted the first version of
    # this analysis: agreement with ELN2017 came out 0.776 when the deduplicated value is 0.822.
    n_rows = len(cl)
    cl = cl[cl["_spec"].str.lower().ne("nan")].drop_duplicates("_spec").reset_index(drop=True)
    cl["_dna"] = cl["dbgap_dnaseq_sample"].astype(str) if "dbgap_dnaseq_sample" in cl.columns else cl["_spec"]
    cl["_karyo"] = cl["karyotype"].map(karyo_flags)

    muts = load_mutations(a.vaf)
    # mutations are keyed by the DNA sample id; map through whichever id matches
    keyed = {}
    for _, r in cl.iterrows():
        for k in (r["_dna"], r["_spec"]):
            if k in muts:
                keyed[r["_spec"]] = muts[k]; break
    hit = len(keyed)

    rows = []
    for _, r in cl.iterrows():
        rr = dict(r); rr["_spec"] = r["_spec"]
        risk, why = classify(rr, keyed, a.vaf)
        rows.append({"specimen": r["_spec"], "ELN2022": risk, "ELN2017": r.get("ELN2017"),
                     "has_mutation_data": r["_spec"] in keyed, "reasons": " | ".join(why)})
    out = pd.DataFrame(rows)
    dst = a.out or os.path.join(ROOT, "labels", "eln2022_beataml_vaf%02d.tsv" % round(100 * a.vaf))
    out.to_csv(dst, sep="\t", index=False)

    print("ELN 2022 inferred at VAF >= %.0f%%" % (100 * a.vaf))
    print("  clinical rows %d -> %d distinct specimens; mutation data for %d of them"
          % (n_rows, len(cl), hit))
    # A specimen with no variant call cannot be assessed for the nine MDS-related genes or for TP53,
    # so it can only ever be under-called adverse. Those are reported separately rather than pooled.
    g = out[out["has_mutation_data"]]
    print("\n  distribution (mutation data available, n=%d): %s"
          % (len(g), g["ELN2022"].value_counts().to_dict()))
    print("  distribution (no mutation data,    n=%d): %s  <- MR/TP53 not assessable"
          % (len(out) - len(g), out[~out["has_mutation_data"]]["ELN2022"].value_counts().to_dict()))
    both = g[g["ELN2017"].isin(["Favorable", "Intermediate", "Adverse"])]
    x = pd.crosstab(both["ELN2017"], both["ELN2022"])
    print("\n  2017 (rows) x 2022 (cols), specimens with mutation data:")
    print(x.to_string())
    agree = float(np.mean(both["ELN2017"].values == both["ELN2022"].values))
    print("\n  agreement with the shipped ELN2017 label: %.3f (n=%d)" % (agree, len(both)))
    print("  wrote %s" % dst)


if __name__ == "__main__":
    main()
