#!/usr/bin/env python3
"""Validate the variant-level (fine) label parser against BeatAML + check non-SNV data + Leucegene cols."""
import os, re, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd

ROOT = r"C:\Users\krog5w\.gemini\antigravity\scratch\AML-multimodal"
BA = os.path.join(ROOT, "data", "external", "beataml")
LG = os.path.join(ROOT, "data", "external", "leucegene")

def fine_labels(sym, hgvsp, vc):
    """protein-change -> fine variant category(ies). Returns list (sample-level biallelic handled later)."""
    hp = str(hgvsp or "").replace("p.", "")
    vc = str(vc or "")
    m = re.match(r"^[A-Za-z*](\d+)", hp)
    pos = int(m.group(1)) if m else None
    out = []
    if sym == "FLT3":
        if pos in (835, 836): out.append("FLT3_TKD_D835/I836")
        elif pos == 676: out.append("FLT3_N676")
        else: out.append("FLT3_other_TKD_or_JM")
    elif sym == "NPM1":
        if "frameshift" in vc or "ins" in vc: out.append("NPM1_exon12_frameshift")
    elif sym == "DNMT3A":
        out.append("DNMT3A_R882" if pos == 882 else "DNMT3A_nonR882")
    elif sym == "IDH1":
        if pos == 132: out.append("IDH1_R132")
    elif sym == "IDH2":
        if pos == 140: out.append("IDH2_R140")
        elif pos == 172: out.append("IDH2_R172")
    elif sym == "NRAS":
        if pos == 12: out.append("NRAS_G12")
        elif pos == 13: out.append("NRAS_G13")
        elif pos == 61: out.append("NRAS_Q61")
    elif sym == "KRAS":
        if pos == 12: out.append("KRAS_G12")
        elif pos == 13: out.append("KRAS_G13")
        elif pos == 61: out.append("KRAS_Q61")
    elif sym == "PTPN11":
        if pos == 76: out.append("PTPN11_E76")
        elif pos == 61: out.append("PTPN11_D61")
        elif pos == 72: out.append("PTPN11_A72")
        else: out.append("PTPN11_other")
    elif sym == "TP53":
        if "missense" in vc and pos and 94 <= pos <= 312: out.append("TP53_hotspot_DBD")
        else: out.append("TP53_LOF/splice/frameshift")
    elif sym == "CEBPA":
        if pos and pos >= 278: out.append("CEBPA_bZIP")
        elif "frameshift" in vc or "stop" in vc: out.append("CEBPA_Nterminal_frameshift/nonsense")
        else: out.append("CEBPA_bZIP")
    elif sym == "RUNX1":
        out.append("RUNX1_LOF")
    elif sym == "SF3B1":
        if pos == 700: out.append("SF3B1_K700")
        elif pos == 666: out.append("SF3B1_K666")
    elif sym == "U2AF1":
        if pos == 34: out.append("U2AF1_S34")
        elif pos in (156, 157): out.append("U2AF1_Q157/R156")
    elif sym == "KIT":
        if pos == 816: out.append("KIT_D816")
        elif pos in (820, 822): out.append("KIT_N822/D820")
    elif sym == "JAK2":
        if pos == 617: out.append("JAK2_V617F")
    elif sym == "WT1":
        out.append("WT1_LOF")
    elif sym == "KMT2D":
        out.append("MLL2_KMT2D_mutation")
    return out

m = pd.read_csv(os.path.join(BA, "mutations.txt"), sep="\t",
                usecols=["dbgap_sample_id", "symbol", "hgvsp_short", "variant_classification"])
# sample-level category -> set of samples
from collections import defaultdict
cat2samp = defaultdict(set)
cebpa_count = defaultdict(int)
for s, g, hp, vc in zip(m["dbgap_sample_id"].astype(str), m["symbol"].astype(str),
                        m["hgvsp_short"], m["variant_classification"]):
    for c in fine_labels(g, hp, vc):
        cat2samp[c].add(s)
    if g == "CEBPA":
        cebpa_count[s] += 1
# CEBPA biallelic/double = >=2 CEBPA variants in a sample
for s, n in cebpa_count.items():
    if n >= 2:
        cat2samp["CEBPA_biallelic_or_double"].add(s)

USER_LIST = ["FLT3_ITD","FLT3_TKD_D835/I836","FLT3_N676","FLT3_other_TKD_or_JM","NPM1_exon12_frameshift",
 "DNMT3A_R882","DNMT3A_nonR882","IDH1_R132","IDH2_R140","IDH2_R172","NRAS_G12","NRAS_G13","NRAS_Q61",
 "KRAS_G12","KRAS_G13","KRAS_Q61","PTPN11_E76","PTPN11_D61","PTPN11_A72","PTPN11_other","CEBPA_bZIP",
 "CEBPA_Nterminal_frameshift/nonsense","CEBPA_biallelic_or_double","TP53_hotspot_DBD","TP53_LOF/splice/frameshift",
 "RUNX1_LOF","SF3B1_K700","SF3B1_K666","U2AF1_S34","U2AF1_Q157/R156","KIT_D816","KIT_N822/D820","JAK2_V617F",
 "WT1_LOF","WT1_rs16754_benign_exclude","KMT2A_PTD","KMT2A_fusion","MLL2_KMT2D_mutation"]

print("==== BeatAML fine-category positive counts (SNV-derived) ====")
for c in USER_LIST:
    n = len(cat2samp.get(c, set()))
    flag = "" if c in ("FLT3_ITD","KMT2A_PTD","KMT2A_fusion","WT1_rs16754_benign_exclude") else ("  <-- <6 (rare)" if n < 6 else "")
    src = " [needs clinical/fusion, not SNV]" if c in ("FLT3_ITD","KMT2A_PTD","KMT2A_fusion") else (" [exclusion tag, not a class]" if "exclude" in c else "")
    print("  %-40s %3d%s%s" % (c, n, src, flag))

print("\n==== BeatAML clinical columns (search FLT3/KMT2A/PTD/NPM1/fusion) ====")
cl = pd.read_excel(os.path.join(BA, "clinical.xlsx"), "summary")
hits = [c for c in cl.columns if re.search(r"FLT3|KMT2A|MLL|PTD|NPM1|fusion|ITD|karyotype|consensus", str(c), re.I)]
for c in hits:
    vc = cl[c].astype(str).value_counts().head(5).to_dict()
    print("  %-28s %s" % (c, {k[:20]: v for k, v in vc.items()}))
print("  [all sheets]:", pd.ExcelFile(os.path.join(BA, "clinical.xlsx")).sheet_names)

print("\n==== files in beataml dir (any fusion/PTD file?) ====")
print(" ", [f for f in os.listdir(BA)])

print("\n==== Leucegene variant columns available ====")
lab = pd.read_excel(os.path.join(LG, "leucegene_labels.xlsx"), "Leucegene Metadata", header=1)
vcols = [c for c in lab.columns if re.search(r"variant|ITD|TDK|TKD|R132|R140|R172|allelic|P95|Q157|S34|G12|hotspot|mutat", str(c), re.I)]
for c in vcols:
    print("   ", c)
