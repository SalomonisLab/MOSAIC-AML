#!/usr/bin/env python3
"""Self-test for pipeline/categories.py. Runs anywhere, no data: python pipeline/_categories_selftest.py"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import categories as cat

# subtypes keep their own rows; only true aliases merge
assert cat.row("FLT3_ITD") != cat.row("FLT3_TKD_D835/I836")
assert cat.row("TP53_hotspot_DBD") != cat.row("TP53_LOF/splice/frameshift")
assert cat.row("inv(16)_CBFB-MYH11") == cat.row("inv16") == "inv16"
assert cat.row("kmt2a") != cat.row("KMT2A-rearrangement"), "two heads, two rows"

# grouping still puts subtypes under their gene, and cytogenetic calls keep their case
assert cat.gene("FLT3-ITD") == cat.gene("FLT3_TKD_D835/I836") == cat.gene("FLT3") == "FLT3"
assert cat.gene("DNMT3A_nonR882") == "DNMT3A"
assert cat.gene("KMT2A-rearrangement") == cat.gene("kmt2a") == cat.gene("KMT2A_fusion") == "KMT2A"
assert [cat.gene(c) for c in ("del7", "del5", "trisomy8", "complex")] == ["del7", "del5", "trisomy8", "complex"]

# an umbrella is left out of a mean only when its subtypes are listed with it
assert cat.excluded("FLT3", ["FLT3", "FLT3-ITD", "FLT3-TKD"]) == "umbrella of FLT3-ITD and FLT3-TKD"
assert cat.excluded("FLT3", ["FLT3", "FLT3-ITD"]) == "umbrella of FLT3-ITD"
assert cat.excluded("FLT3", ["FLT3", "NPM1"]) is None
assert cat.excluded("KMT2A-rearrangement", ["kmt2a", "KMT2A-rearrangement"]).startswith("second head")
assert cat.excluded("KMT2A-rearrangement", ["KMT2A-rearrangement"]) is None
assert cat.excluded("kmt2a", ["kmt2a", "KMT2A-rearrangement"]) is None

# bulk-to-multimodal matches: only the lesions the domain review called the same or nearly so
assert cat.match("FLT3_ITD") == ("FLT3-ITD", "same")
assert cat.match("NPM1_exon12_frameshift") == ("NPM1", "near")
assert cat.match("KIT_D816") == (None, None), "KIT holds more than D816"
assert cat.match("PTPN11_other") == (None, None), "PTPN11_other excludes the hotspots"
assert cat.match("FLT3_other_TKD_or_JM") == (None, None)

# Wilson 95% interval, against values computed by hand
lo, hi = cat.wilson(0, 5)
assert lo == 0.0 and abs(hi - 0.4345) < 1e-3, (lo, hi)
lo, hi = cat.wilson(112, 158)
assert abs(lo - 0.6333) < 1e-3 and abs(hi - 0.7740) < 1e-3, (lo, hi)
assert all(x != x for x in cat.wilson(0, 0))

# AUROC: perfect, inverted, ties count one half, undefined without both classes
assert cat.auroc([0.9, 0.8], [0.1, 0.2]) == 1.0
assert cat.auroc([0.1], [0.9]) == 0.0
assert cat.auroc([0.5, 0.9], [0.5]) == 0.75
assert cat.auroc([], [0.3]) is None and cat.auroc([0.3], []) is None

print("CATEGORIES SELFTEST OK")
