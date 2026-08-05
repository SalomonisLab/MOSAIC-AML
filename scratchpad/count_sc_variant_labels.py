#!/usr/bin/env python3
"""Parse labels/vaf_per_sample.tsv HGVS -> variant categories; count sc-cohort positives per category
among genotyped samples (index of mutation_matrix_explicit_v2.tsv). Decides which clear the >=8 floor."""
import os, re, sys
import pandas as pd
from collections import defaultdict
ROOT = r"C:\Users\krog5w\.gemini\antigravity\scratch\AML-multimodal"
sys.path.insert(0, os.path.join(ROOT, "pipeline"))
import bulk_external as BE   # reuse the _fine parser + taxonomy

# authoritative genotyped samples + per-gene present/absent
MM = pd.read_csv(os.path.join(ROOT, "pipeline", "mutation_matrix_explicit_v2.tsv"), sep="\t", index_col=0)
genotyped = set(MM.index)
gene_present = {c[4:]: set(MM.index[MM[c] > 0]) for c in MM.columns if c.startswith("mut_")}
print("genotyped samples: %d | genes in explicit matrix: %d" % (len(genotyped), len(gene_present)))

vaf = pd.read_csv(os.path.join(ROOT, "labels", "vaf_per_sample.tsv"), sep="\t")
vaf = vaf[vaf["in_cohort"] == 1]

def parse_protein(hgvs, gene):
    """hgvs like 'ASXL1.L815P 100%' / 'ASXL1: G646Wfs 39%' -> protein token 'L815P'."""
    s = str(hgvs)
    # drop the gene prefix + separators
    s = re.sub(r"^%s\s*[.: ]+" % re.escape(gene), "", s)
    s = s.split()[0] if s.split() else s          # first whitespace token = protein change
    return s

def pseudo_vc(prot):
    p = prot.lower()
    if "fs" in p: return "frameshift_variant"
    if "*" in p or "ter" in p: return "stop_gained"
    if "del" in p: return "inframe_deletion"
    if "ins" in p or "dup" in p or "itd" in p: return "inframe_insertion"
    if "splice" in p: return "splice_acceptor_variant"
    return "missense_variant"

cat_samp = defaultdict(set)
cebpa_ct = defaultdict(int)
flt3_itd = gene_present.get("FLT3-ITD", set())     # already an explicit split column if present
unparsed = defaultdict(int)
for _, r in vaf.iterrows():
    g = str(r["gene"]); s = str(r["sample_key"]); prot = parse_protein(r["hgvs"], g)
    if s not in genotyped:
        continue
    if g in BE.SPLIT_GENES:
        vc = pseudo_vc(prot)
        for c in BE._fine(g, "p." + prot, vc):
            cat_samp[c].add(s)
        if g == "CEBPA":
            cebpa_ct[s] += 1
        if not BE._fine(g, "p." + prot, vc):
            unparsed[g] += 1
for s, n in cebpa_ct.items():
    if n >= 2: cat_samp["CEBPA_biallelic_or_double"].add(s)
# FLT3_ITD from explicit column
cat_samp["FLT3_ITD"] |= flt3_itd

print("\n== variant categories present in sc cohort (>=1 positive), sorted by count ==")
rows = []
for c, ss in cat_samp.items():
    rows.append((len(ss & genotyped), c))
for n, c in sorted(rows, reverse=True):
    flag = "" if n >= 8 else ("  (<8: not trainable, needs 8)" if n >= 1 else "")
    print("  %-38s %3d%s" % (c, n, flag))

print("\n== gene-level (non-split) genes with >=8 present ==")
for g in BE.GENE_DRIVERS:
    n = len(gene_present.get(g, set()))
    if n >= 1:
        print("  %-12s %3d%s" % (g, n, "" if n >= 8 else "  (<8)"))
if unparsed:
    print("\n(unparsed split-gene variants:", dict(unparsed), ")")
