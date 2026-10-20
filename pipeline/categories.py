"""Mutation-category names, shared by every script that reports results by category.

Four scripts each carried their own canon(), which cut a name at the first
"_" or "-" and so merged distinct subtypes (FLT3 ITD and TKD, TP53 hotspot and
LOF) into one key, kept whichever came last, and counted aliases twice. The
rules here are explicit instead, each with the reason it holds.

  gene(name)        the gene a category belongs to, for grouping rows
  row(name)         the row a category is reported under (aliases merged)
  excluded(name, names)
                    why a category is listed but left out of a mean, or None
  match(bulk_name)  the multimodal row naming the same lesion, and how closely
"""
import math

# Two names for one call: reported once, under the first name. inv16 and
# inv(16)_CBFB-MYH11 agree on every held-out patient.
ALIASES = {"inv(16)_CBFB-MYH11": "inv16"}

# A category whose own subtypes are also reported. It is listed, but a mean
# that included it would count those patients twice.
UMBRELLAS = {"FLT3": ("FLT3-ITD", "FLT3-TKD")}

# A second model head for a lesion another head reports. kmt2a is the head the
# arbiter uses for the KMT2A call (amlmm/arbiter.py ANCHOR_PRIORITY); the two
# disagree on two held-out patients, so both are listed and one is averaged.
SECOND_HEADS = {"KMT2A-rearrangement": "kmt2a"}

# Bulk (BeatAML-trained) category -> (multimodal category, closeness).
# "same": one lesion by definition or by the bulk labeling code.
# "near": the bulk subtype is nearly all of the gene-level lesion in AML.
# Not matched, because the multimodal category clearly holds more: KIT vs
# KIT_D816 (about 41% of KIT mutations in CBF-AML fall outside D816) and
# PTPN11 vs PTPN11_other (the bulk category excludes the hotspots).
MATCH = {
    "FLT3_ITD": ("FLT3-ITD", "same"),
    "KMT2A_fusion": ("kmt2a", "same"),
    "RUNX1_LOF": ("RUNX1", "same"),      # bulk_external.py takes every RUNX1 variant
    "WT1_LOF": ("WT1", "same"),          # likewise every WT1 variant
    "NPM1_exon12_frameshift": ("NPM1", "near"),   # about 1.3% of NPM1-mutant AML is outside exon 12
    "IDH1_R132": ("IDH1", "near"),       # about 4% of IDH1 mutations are not R132
    "FLT3_TKD_D835/I836": ("FLT3-TKD", "near"),   # D835/I836 is the clinical TKD assay
}


def gene(name):
    """The gene a category belongs to, for grouping and ordering rows."""
    s = str(name)
    low = s.lower()
    if "inv16" in low or "inv(16)" in low or "cbfb" in low:
        return "inv16"
    if "kmt2a" in low:
        return "KMT2A"
    if low in ("del5", "del7", "trisomy8", "complex"):
        return low
    return s.split("_")[0].split("-")[0].upper()


def row(name):
    """The row a category is reported under."""
    return ALIASES.get(str(name), str(name))


def excluded(name, names):
    """Why `name` is listed but kept out of a mean over `names`, or None."""
    names = set(names)
    subs = [s for s in UMBRELLAS.get(name, ()) if s in names]
    if subs:
        return "umbrella of %s" % " and ".join(subs)
    other = SECOND_HEADS.get(name)
    if other in names:
        return "second head; the mean uses %s" % other
    return None


def match(bulk_name):
    """(multimodal row, 'same' | 'near') for a bulk category, else (None, None).

    A bulk category with exactly a multimodal category's name matches it.
    """
    if bulk_name in MATCH:
        return MATCH[bulk_name]
    return None, None


def identical_name(bulk_name, multimodal_names):
    return bulk_name if bulk_name in set(multimodal_names) else None


def auroc(pos, neg):
    """Area under the ROC curve from scores of positives and negatives
    (Mann-Whitney, ties count one half); None without both classes."""
    if not pos or not neg:
        return None
    wins = 0.0
    for p in pos:
        for n in neg:
            wins += 1.0 if p > n else 0.5 if p == n else 0.0
    return wins / (len(pos) * len(neg))


def wilson(k, n, z=1.959964):
    """Wilson score 95% interval for k successes in n trials; (nan, nan) if n is 0."""
    if not n:
        return float("nan"), float("nan")
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0.0, c - h), min(1.0, c + h)
