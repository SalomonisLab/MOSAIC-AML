#!/usr/bin/env python3
"""Harmonize three cohorts into ONE gene-expression space + common driver labels, for the
bulk / bulk+sc mutation-prediction bake-off (BeatAML2, Leucegene, and our scRNA-as-bulk).

  - BeatAML2 : norm_exp (22.8k genes x ~700 RNA samples) + WES mutations + curated clinical + VAF.
  - Leucegene: AltAnalyze steady-state (19k genes x 367 samples) + dbGAP variant labels (KNOWN/0/UNK).
  - our sc   : whole-sample bulk-equivalent (sum all cells -> CP10k+log1p) from the atlas + our labels.

Labels are VARIANT-LEVEL (per Nathan's spec): genes with distinct functional hotspots are split into
sub-categories (FLT3_ITD vs TKD vs N676; DNMT3A_R882 vs nonR882; NRAS_G12/G13/Q61; TP53 hotspot-DBD vs
LOF; etc.) parsed from BeatAML hgvsp_short + variant_classification + clinical fusions. Genes without a
meaningful functional split stay gene-level. BeatAML CV-OOF is the primary metric (covers every category
with >=6 positives); Leucegene / single-cell externally validate whichever categories they carry calls for.

Everything keyed on Ensembl gene IDs (ENSG). Per-dataset z-scoring handles the cross-cohort batch diff.

Run:  python bulk_external.py     # prints the harmonization feasibility report
"""
import os, sys, re, warnings
from collections import defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
BA = os.path.join(ROOT, "data", "external", "beataml")
LG = os.path.join(ROOT, "data", "external", "leucegene")

# --- label taxonomy -------------------------------------------------------------------------------
# genes split into variant-level categories (distinct biology); base gene-level label is NOT used for these.
SPLIT_GENES = {"FLT3", "NPM1", "DNMT3A", "IDH1", "IDH2", "NRAS", "KRAS", "PTPN11", "CEBPA", "TP53",
               "RUNX1", "SF3B1", "U2AF1", "KIT", "JAK2", "WT1", "KMT2D"}
# fine categories built from SNV protein change (hgvsp_short) + variant_classification
FINE_SNV = ["FLT3_TKD_D835/I836", "FLT3_N676", "FLT3_other_TKD_or_JM", "NPM1_exon12_frameshift",
            "DNMT3A_R882", "DNMT3A_nonR882", "IDH1_R132", "IDH2_R140", "IDH2_R172",
            "NRAS_G12", "NRAS_G13", "NRAS_Q61", "KRAS_G12", "KRAS_G13", "KRAS_Q61",
            "PTPN11_E76", "PTPN11_D61", "PTPN11_A72", "PTPN11_other", "CEBPA_bZIP",
            "CEBPA_Nterminal_frameshift/nonsense", "CEBPA_biallelic_or_double", "TP53_hotspot_DBD",
            "TP53_LOF/splice/frameshift", "RUNX1_LOF", "SF3B1_K700", "SF3B1_K666", "U2AF1_S34",
            "U2AF1_Q157/R156", "KIT_D816", "KIT_N822/D820", "JAK2_V617F", "WT1_LOF", "MLL2_KMT2D_mutation"]
# categories from curated clinical / fusion calls (not SNV)
FINE_CLIN = ["FLT3_ITD", "KMT2A_fusion", "KMT2A_PTD"]
# genes kept whole (no functional split): gene-level presence/absence
GENE_DRIVERS = ["TET2", "ASXL1", "SRSF2", "STAG2", "BCOR", "GATA2", "PHF6", "EZH2", "NF1", "RAD21",
                "CBL", "BCORL1", "SMC3", "SMC1A", "ZRSR2", "IKZF1", "CSF3R", "SUZ12", "SETBP1",
                "ETV6", "CREBBP"]
DRIVERS = FINE_SNV + FINE_CLIN + GENE_DRIVERS
_SNV_CATS = set(FINE_SNV) | set(GENE_DRIVERS)   # categories that take "0" for any WES'd sample


def _fine(sym, hgvsp, vc):
    """Protein change -> fine variant category(ies). BeatAML hgvsp_short like 'p.D835Y', 'p.R882H'."""
    hp = str(hgvsp or "").replace("p.", "")
    vc = str(vc or "")
    m = re.match(r"^[A-Za-z*](\d+)", hp)
    pos = int(m.group(1)) if m else None
    o = []
    if sym == "FLT3":
        if pos in (835, 836): o.append("FLT3_TKD_D835/I836")            # activation-loop TKD
        elif pos == 676: o.append("FLT3_N676")
        else: o.append("FLT3_other_TKD_or_JM")                          # ITD comes from clinical, not here
    elif sym == "NPM1":
        if "frameshift" in vc or "ins" in vc: o.append("NPM1_exon12_frameshift")
    elif sym == "DNMT3A":
        o.append("DNMT3A_R882" if pos == 882 else "DNMT3A_nonR882")
    elif sym == "IDH1":
        if pos == 132: o.append("IDH1_R132")
    elif sym == "IDH2":
        if pos == 140: o.append("IDH2_R140")
        elif pos == 172: o.append("IDH2_R172")
    elif sym == "NRAS":
        if pos == 12: o.append("NRAS_G12")
        elif pos == 13: o.append("NRAS_G13")
        elif pos == 61: o.append("NRAS_Q61")
    elif sym == "KRAS":
        if pos == 12: o.append("KRAS_G12")
        elif pos == 13: o.append("KRAS_G13")
        elif pos == 61: o.append("KRAS_Q61")
    elif sym == "PTPN11":
        if pos == 76: o.append("PTPN11_E76")
        elif pos == 61: o.append("PTPN11_D61")
        elif pos == 72: o.append("PTPN11_A72")
        else: o.append("PTPN11_other")
    elif sym == "TP53":
        if "missense" in vc and pos and 94 <= pos <= 312: o.append("TP53_hotspot_DBD")   # DNA-binding domain
        else: o.append("TP53_LOF/splice/frameshift")
    elif sym == "CEBPA":
        if pos and pos >= 278: o.append("CEBPA_bZIP")                   # C-terminal bZIP domain
        elif "frameshift" in vc or "stop" in vc: o.append("CEBPA_Nterminal_frameshift/nonsense")
        else: o.append("CEBPA_bZIP")
    elif sym == "RUNX1":
        o.append("RUNX1_LOF")
    elif sym == "SF3B1":
        if pos == 700: o.append("SF3B1_K700")
        elif pos == 666: o.append("SF3B1_K666")
    elif sym == "U2AF1":
        if pos == 34: o.append("U2AF1_S34")
        elif pos in (156, 157): o.append("U2AF1_Q157/R156")
    elif sym == "KIT":
        if pos == 816: o.append("KIT_D816")
        elif pos in (820, 822): o.append("KIT_N822/D820")
    elif sym == "JAK2":
        if pos == 617: o.append("JAK2_V617F")
    elif sym == "WT1":
        o.append("WT1_LOF")                                            # benign rs16754 already excluded by curated calls
    elif sym == "KMT2D":
        o.append("MLL2_KMT2D_mutation")
    return o


def load_beataml():
    ex = pd.read_csv(os.path.join(BA, "norm_exp.txt"), sep="\t")
    ens = ex["stable_id"].astype(str).values
    sym2ens = {}
    for e, s in zip(ens, ex["display_label"].astype(str)):
        sym2ens.setdefault(s, e)                                       # first ENSG per symbol
    scols = list(ex.columns[4:])                                       # BA####R RNA samples
    X = ex[scols].T.copy(); X.columns = ens; X.index = scols
    X = X.loc[:, ~pd.Index(ens).duplicated()]
    X = np.power(2.0, X.astype(float))                                 # BeatAML norm_exp is log2 -> linearize

    mut = pd.read_csv(os.path.join(BA, "mutations.txt"), sep="\t",
                      usecols=["dbgap_sample_id", "symbol", "hgvsp_short", "variant_classification"])
    cl = pd.read_excel(os.path.join(BA, "clinical.xlsx"), "summary")
    d2r = {str(d): str(r) for d, r in zip(cl["dbgap_dnaseq_sample"], cl["dbgap_rnaseq_sample"]) if pd.notna(r)}

    lab = pd.DataFrame(index=X.index, columns=DRIVERS, dtype=float)
    # any sample with WES (DNA mapped to an RNA sample we have) -> 0 by default for SNV-derived categories
    wes_rna = set(d2r.get(str(d)) for d in mut["dbgap_sample_id"].unique()) & set(X.index)
    for c in _SNV_CATS:
        lab.loc[list(wes_rna), c] = 0.0
    cebpa = defaultdict(int)
    for s, g, hp, vc in zip(mut["dbgap_sample_id"].astype(str), mut["symbol"].astype(str),
                            mut["hgvsp_short"], mut["variant_classification"]):
        r = d2r.get(s)
        if r not in lab.index:
            continue
        if g in GENE_DRIVERS:
            lab.loc[r, g] = 1.0
        for c in _fine(g, hp, vc):
            if c in lab.columns:
                lab.loc[r, c] = 1.0
        if g == "CEBPA":
            cebpa[r] += 1
    for r, n in cebpa.items():                                         # biallelic/double = >=2 CEBPA variants
        if r in lab.index and n >= 2:
            lab.loc[r, "CEBPA_biallelic_or_double"] = 1.0
    # curated clinical: FLT3-ITD (219 pos) + KMT2A fusions (~53). KMT2A_PTD: no BeatAML source -> NaN.
    for r, v in zip(cl["dbgap_rnaseq_sample"].astype(str), cl["FLT3-ITD"].astype(str)):
        if r in lab.index:
            lab.loc[r, "FLT3_ITD"] = 1.0 if v == "positive" else (0.0 if v == "negative" else np.nan)
    for r, v in zip(cl["dbgap_rnaseq_sample"].astype(str), cl["consensusAMLFusions"].astype(str)):
        if r in lab.index:
            lab.loc[r, "KMT2A_fusion"] = 1.0 if "KMT2A" in v.upper() else 0.0
    return X, lab, sym2ens


def _code(v):
    """Leucegene label cell -> 1/0/NaN. Handles KNOWN/0/UNK and 1/0/TRUE/FALSE/WT and residue strings."""
    s = str(v).strip().upper()
    if s in ("KNOWN", "1", "1.0", "TRUE", "YES", "MUT", "POS", "POSITIVE"):
        return 1.0
    if s in ("0", "0.0", "FALSE", "NO", "WT", "NEG", "NEGATIVE", "ABSENT", "-"):
        return 0.0
    if s in ("", "NAN", "NA", "UNK", "UNKNOWN", "NONE", "."):
        return np.nan
    return 1.0                                                         # a variant descriptor present -> positive


def load_leucegene():
    fp = os.path.join(LG, "exp.Leucegene-steady-state.txt")
    ex = pd.read_csv(fp, sep="\t")
    gcol = ex.columns[0]
    ens = ex[gcol].astype(str).values
    scols = list(ex.columns[1:])
    samp = [c.replace(".bed", "").split("_")[0] for c in scols]
    X = ex[scols].T.copy(); X.columns = ens; X.index = samp
    X = X.loc[:, ~pd.Index(ens).duplicated()]; X = X.loc[~pd.Index(samp).duplicated()]

    lab = pd.read_excel(os.path.join(LG, "leucegene_labels.xlsx"), "Leucegene Metadata", header=1)
    lab = lab.dropna(subset=["SRA_Sample_ID"])
    lab.index = lab["SRA_Sample_ID"].astype(str).str.split("_").str[0]
    # map each label category to the Leucegene column that carries residue-level (or gene-level) calls
    COL = {
        # fine categories Leucegene resolves at residue/event level:
        "FLT3_ITD": "FLT3-ITD_variant", "FLT3_TKD_D835/I836": "FLT3-TDK_variant",
        "FLT3_other_TKD_or_JM": "Mutation.FLT3-Other", "IDH1_R132": "IDH1-R132_variants",
        "IDH2_R140": "IDH2-R140Q_variants", "U2AF1_S34": "U2AF1-S34_variants",
        "U2AF1_Q157/R156": "U2AF1-Q157_variants", "NPM1_exon12_frameshift": "NPM1_variant",
        "CEBPA_biallelic_or_double": "CEBPA-biallelic_variants", "KMT2A_PTD": "MLL-PTD_variant",
        # split-gene fine cats where Leucegene only has a gene-level call -> approximate to the dominant hotspot:
        "JAK2_V617F": "Mutation.JAK2", "RUNX1_LOF": "Mutation.RUNX1", "WT1_LOF": "WT1_variants",
        # gene-level (non-split) drivers:
        "TET2": "TET2_variants", "ZRSR2": "ZRSR2_variants", "SETBP1": "Mutation.SETBP1",
        "RAD21": "Mutation.RAD21", "EZH2": "Mutation.EZH2", "PHF6": "Mutation.PHF6",
        "CBL": "Mutation.CBL", "ETV6": "Mutation.ETV6", "GATA2": "Mutation.GATA2",
        "IKZF1": "Mutation.IKZF1", "CSF3R": "Mutation.CSF3R", "STAG2": "Mutation.STAG2",
        "SMC1A": "Mutation.SMC1A", "SMC3": "Mutation.SMC3", "BCOR": "Mutation.BCOR",
        "ASXL1": "Mutation.ASXL1", "SRSF2": "SRSF2-P95_variants",
    }
    L = pd.DataFrame(index=lab.index, columns=DRIVERS, dtype=float)
    for drv in DRIVERS:
        if drv in COL and COL[drv] in lab.columns:
            L[drv] = lab[COL[drv]].map(_code)
    L = L[~L.index.duplicated()].reindex(X.index)                      # align labels to expression samples
    return X, L


def load_sc_bulk(sym2ens):
    import bulk_features as BF
    from amlmm.context import build_context, Config
    from amlmm import discovery as D, genetics
    ctx = build_context(Config(run_id="bulk_external"))
    B = BF.bulk_rna_matrix(ctx)                                         # samples x gene symbols
    cols = [c for c in B.columns if c in sym2ens]
    X = B[cols].copy(); X.columns = [sym2ens[c] for c in cols]         # -> ENSG
    X = X.loc[:, ~pd.Index(X.columns).duplicated()]
    X = np.expm1(X)                                                     # our sc bulk is CP10k+log1p -> linearize
    M = ctx.tables.get("mutations") or genetics.build_mutation_matrix(ctx)
    _m01 = {"present": 1.0, "absent": 0.0}
    L = pd.DataFrame(index=X.index, columns=DRIVERS, dtype=float)
    # our single-cell mutation matrix is gene-level (+ FLT3-ITD split); map what we have, rest stay NaN
    sc_name = {"FLT3_ITD": "mut_FLT3-ITD"}                              # gene-level non-split -> mut_<gene>
    for g in GENE_DRIVERS:
        sc_name[g] = "mut_" + g
    for drv in DRIVERS:
        flag = sc_name.get(drv)
        if flag and flag in M.columns:
            L[drv] = D._labels_for_field_raw(ctx, flag).reindex(X.index).map(_m01)
    return X, L, set(ctx.holdout), ctx


def _clog_z(Xlin):
    clog = np.log2(Xlin.clip(lower=0) + 1.0)                            # common log2 scale for all cohorts
    mu, sd = clog.mean(0), clog.std(0).replace(0, 1.0)
    return (clog - mu) / sd                                             # per-dataset z-scored (batch handling)


def harmonize():
    baX, baL, sym2ens = load_beataml()                                 # all three returned LINEAR
    lgX, lgL = load_leucegene()
    scX, scL, hold, ctx = load_sc_bulk(sym2ens)
    common = sorted(set(baX.columns) & set(lgX.columns) & set(scX.columns))
    out = {"genes": common, "holdout": hold, "ctx": ctx}
    for name, X, L in [("beataml", baX, baL), ("leucegene", lgX, lgL), ("sc", scX, scL)]:
        xl = X[common]
        out[name] = {"Xlin": xl, "Xz": _clog_z(xl), "L": L}            # Xlin = linear (fold FS); Xz = classifier input
    return out


if __name__ == "__main__":
    H = harmonize()
    g = H["genes"]
    print("common ENSG genes:", len(g), "| categories:", len(DRIVERS))
    print("%-38s %-14s %-14s %-14s" % ("category", "beataml", "leucegene", "sc"))
    for d in DRIVERS:
        cells = []
        for c in ["beataml", "leucegene", "sc"]:
            L = H[c]["L"]
            cells.append("%3d/%4d" % (int((L[d] == 1).sum()), int(L[d].notna().sum())))
        star = "" if int((H["beataml"]["L"][d] == 1).sum()) >= 6 else "  (BA<6)"
        print("%-38s %-14s %-14s %-14s%s" % (d, cells[0], cells[1], cells[2], star))
    print("(cells = positives / labeled ; sc holdout n =", len(H["holdout"] & set(H["sc"]["Xz"].index)), ")")
