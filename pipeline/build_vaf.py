#!/usr/bin/env python3
"""Build the per-mutation VAF (variant allele frequency) resource for MOSAIC-AML.

VAF = the fraction of sequencing reads carrying the variant ~ clonal burden of the mutation. It is NOT
produced by our single-cell pipeline (which yields binary present/absent labels); it comes from external
DNA panels. This script harvests whatever VAF we have from the banked cohort files, maps each value to a
deployed driver, and writes:

  labels/vaf_per_sample.tsv   -- canonical long table: one row per (sample, variant) with VAF + provenance
  gui/vaf_by_mutation.json    -- per-DEPLOYED-mutation aggregate, with an explicit slot for EVERY mutation:
                                   status = has_vaf | awaiting (SNV, no data yet) | not_applicable (cytogenetic)

To add a cohort: drop its VAF file in labels/vaf_sources/, extend load_* below, and re-run. Mutations whose
slot is "awaiting" fill in automatically as soon as matching VAF appears — no schema change needed.

Run:  python build_vaf.py
"""
import os, sys, json, re, time, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
warnings.filterwarnings("ignore")
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SRC = os.path.join(ROOT, "labels", "vaf_sources")
CARD = os.path.join(HERE, "model_card.json")

# cytogenetic / structural drivers — VAF is not defined for these (karyotype/FISH events, not point mutations)
CYTO = {"complex", "del5", "del7", "inv16", "trisomy8", "kmt2a",
        "KMT2A-rearrangement", "inv(16)_CBFB-MYH11"}
# SNV drivers whose VAF gene == the driver name (FLT3 handled specially for ITD/TKD split)
GENES_OF_INTEREST = ["FLT3", "NPM1", "IDH1", "IDH2", "NRAS", "KRAS", "TET2", "DNMT3A", "TP53",
                     "RUNX1", "WT1", "PTPN11", "KIT", "ASXL1", "CEBPA", "GATA2"]


def flt3_subtype(text):
    t = str(text).upper()
    if "ITD" in t or "DUP" in t or ("INS" in t and "835" not in t):
        return "FLT3-ITD"
    if "835" in t or "836" in t or "TKD" in t or "D835" in t:
        return "FLT3-TKD"
    return None


def load_meta_cchmc(rows):
    fp = os.path.join(SRC, "Meta-CCHMC.xlsx")
    if not os.path.exists(fp):
        return
    for sheet, tp, vcol in [("Pretreatment samples", "pretreatment", "VAF"),
                            ("Relapse samples", "relapse", "VAF (at relapse)")]:
        try:
            df = pd.read_excel(fp, sheet)
        except Exception:
            continue
        for _, r in df.iterrows():
            g = str(r.get("Gene", "")).strip()
            v = r.get(vcol)
            if g not in GENES_OF_INTEREST or pd.isna(v):
                continue
            try:
                v = float(v)
            except Exception:
                continue
            if not (0 < v <= 1):
                continue
            sid = str(r.get("scRNA-seq ID", "")).strip()
            hg = str(r.get("Mutation", "")).strip()
            drivers = [g]
            if g == "FLT3":
                st = flt3_subtype(hg)
                if st:
                    drivers.append(st)
            for d in drivers:
                rows.append({"sample_key": ("CCHMC::" + sid) if sid else "", "driver": d, "gene": g,
                             "vaf": round(v, 4), "hgvs": hg[:80], "timepoint": tp, "source": "Meta-CCHMC"})


def load_variant_detail(rows):
    fp = os.path.join(SRC, "AML_harmonized_metadata.xlsx")
    if not os.path.exists(fp):
        return
    try:
        df = pd.read_excel(fp, "Variant_Detail")
    except Exception:
        return
    genepat = re.compile(r"(" + "|".join(GENES_OF_INTEREST) + r")", re.I)
    pctpat = re.compile(r"(\d{1,3}(?:\.\d+)?)\s*%")
    for _, r in df.iterrows():
        ds, sm = str(r.get("Dataset", "")).strip(), str(r.get("Sample", "")).strip()
        key = "%s::%s" % (ds, sm)
        txt = str(r.get("Variant (HGVS/desc)", ""))
        for seg in re.split(r"[,;]", txt):
            gm = genepat.search(seg)
            pm = pctpat.search(seg)
            if not (gm and pm):
                continue
            g = gm.group(1).upper()
            if g not in GENES_OF_INTEREST:
                continue
            v = float(pm.group(1)) / 100.0
            if not (0 < v <= 1):
                continue
            drivers = [g]
            if g == "FLT3":
                st = flt3_subtype(seg)
                if st:
                    drivers.append(st)
            for d in drivers:
                rows.append({"sample_key": key, "driver": d, "gene": g, "vaf": round(v, 4),
                             "hgvs": seg.strip()[:80], "timepoint": "diagnosis", "source": "Variant_Detail"})


def quant(vals, q):
    if not vals:
        return None
    s = sorted(vals); i = q * (len(s) - 1); lo = int(i); frac = i - lo
    return round(s[lo] + (s[min(lo + 1, len(s) - 1)] - s[lo]) * frac, 3)


def main():
    muts = list(json.load(open(CARD))["mutations"].keys()) if os.path.exists(CARD) else []
    rows = []
    load_meta_cchmc(rows)
    load_variant_detail(rows)
    df = pd.DataFrame(rows)
    # our cohort sample keys (flag whether a VAF sample is in the single-cell atlas)
    try:
        from amlmm.context import build_context, Config
        ours = set(build_context(Config(run_id="vafbuild")).tables["samples"].index)
    except Exception:
        ours = set()
    if not df.empty:
        df["in_cohort"] = df["sample_key"].isin(ours).astype(int)
        df = df.sort_values(["driver", "sample_key"]).drop_duplicates()
        df.to_csv(os.path.join(ROOT, "labels", "vaf_per_sample.tsv"), sep="\t", index=False)

    # per-mutation aggregate: one representative VAF per (sample, driver) = max (dominant clone)
    out = {"meta": {"built": time.strftime("%Y-%m-%d %H:%M"),
                    "definition": "VAF = variant allele fraction (~clonal burden). External DNA-panel data, "
                                  "not produced by the single-cell pipeline.",
                    "sources": ["Meta-CCHMC.xlsx (structured)", "AML_harmonized_metadata.xlsx / Variant_Detail (parsed)"],
                    "status_key": "has_vaf = VAF available | awaiting = SNV driver, no VAF yet (slot ready) | "
                                  "not_applicable = cytogenetic/structural lesion, VAF undefined"},
           "mutations": {}}
    n_has = n_wait = n_na = 0
    for m in muts:
        if m in CYTO:
            out["mutations"][m] = {"kind": "cytogenetic", "status": "not_applicable", "median": None,
                                   "n_samples": 0, "note": "structural/karyotype lesion — VAF not defined"}
            n_na += 1
            continue
        sub = df[df["driver"] == m] if not df.empty else pd.DataFrame()
        if len(sub):
            # per-sample dominant VAF
            per = sub.groupby("sample_key")["vaf"].max() if "sample_key" in sub else sub["vaf"]
            vals = [float(v) for v in per.values]
            bycoh = {}
            for k in sub["sample_key"]:
                c = str(k).split("::")[0] if "::" in str(k) else "?"
                bycoh[c] = bycoh.get(c, 0) + 1
            out["mutations"][m] = {"kind": "snv", "status": "has_vaf", "n_samples": len(vals),
                                   "median": quant(vals, 0.5), "q1": quant(vals, 0.25), "q3": quant(vals, 0.75),
                                   "min": round(min(vals), 3), "max": round(max(vals), 3),
                                   "values": [round(v, 3) for v in sorted(vals)],
                                   "by_cohort": bycoh}
            n_has += 1
        else:
            out["mutations"][m] = {"kind": "snv", "status": "awaiting", "median": None, "n_samples": 0,
                                   "note": "awaiting VAF — slot ready; add a cohort file to labels/vaf_sources/ and re-run"}
            n_wait += 1
    out["meta"].update({"n_has_vaf": n_has, "n_awaiting": n_wait, "n_cytogenetic": n_na,
                        "n_variant_rows": int(len(df))})
    json.dump(out, open(os.path.join(ROOT, "gui", "vaf_by_mutation.json"), "w"), indent=1)
    print("wrote labels/vaf_per_sample.tsv (%d variant rows) and gui/vaf_by_mutation.json" % len(df))
    print("per-mutation VAF slots: %d has_vaf, %d awaiting, %d cytogenetic(N/A)" % (n_has, n_wait, n_na))
    for m in muts:
        e = out["mutations"][m]
        med = ("median %.2f n=%d" % (e["median"], e["n_samples"])) if e["status"] == "has_vaf" else e["status"]
        print("  %-22s %s" % (m, med))


if __name__ == "__main__":
    main()
