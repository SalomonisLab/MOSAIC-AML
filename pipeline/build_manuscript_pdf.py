#!/usr/bin/env python3
"""Typeset the MOSAIC-AML manuscript to PDF with the figures and tables embedded.

No LaTeX/pandoc on this host, so the document is composed directly with ReportLab Platypus:
markdown headings, paragraphs, block quotes, tables and inline emphasis are mapped to flowables, the
five journal-format figures are placed at their callouts with legends, and the tables are rendered from
the generated TSVs so they cannot drift from the analyses.

  python build_manuscript_pdf.py   -> deliverables/MOSAIC-AML_manuscript.pdf
"""
import os, re, csv, glob
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import mm
from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.enums import TA_JUSTIFY, TA_CENTER
from reportlab.platypus import (BaseDocTemplate, PageTemplate, Frame, Paragraph, Spacer,
                                Image, Table, TableStyle, KeepTogether, PageBreak)
from PIL import Image as PILImage

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D    = os.path.join(ROOT, "deliverables")
FIGD = os.path.join(D, "figures"); TABD = os.path.join(D, "tables")
SRC  = os.path.join(D, "MOSAIC-AML_manuscript.md")
OUT  = os.path.join(D, "MOSAIC-AML_manuscript.pdf")

INK = colors.HexColor("#1a1a19"); MUT = colors.HexColor("#5b5b54")
ACC = colors.HexColor("#12507a"); LINE = colors.HexColor("#c9c9c2")
BODY_W = A4[0] - 40*mm

S = {
 "title": ParagraphStyle("title", fontName="Helvetica-Bold", fontSize=15, leading=19,
                         textColor=INK, spaceAfter=5),
 "sub":   ParagraphStyle("sub", fontName="Helvetica", fontSize=9, leading=12, textColor=MUT, spaceAfter=2),
 "h1":    ParagraphStyle("h1", fontName="Helvetica-Bold", fontSize=12.5, leading=15, textColor=ACC,
                         spaceBefore=11, spaceAfter=5),
 "h2":    ParagraphStyle("h2", fontName="Helvetica-Bold", fontSize=10, leading=13, textColor=INK,
                         spaceBefore=8, spaceAfter=3),
 "body":  ParagraphStyle("body", fontName="Helvetica", fontSize=8.9, leading=12.4, textColor=INK,
                         alignment=TA_JUSTIFY, spaceAfter=5),
 "quote": ParagraphStyle("quote", fontName="Helvetica-Oblique", fontSize=8.9, leading=12.4,
                         textColor=INK, leftIndent=10, rightIndent=10, spaceBefore=3, spaceAfter=6,
                         borderPadding=4, backColor=colors.HexColor("#f4f4ef")),
 "cap":   ParagraphStyle("cap", fontName="Helvetica", fontSize=7.6, leading=10, textColor=MUT,
                         alignment=TA_JUSTIFY, spaceBefore=3, spaceAfter=9),
 "li":    ParagraphStyle("li", fontName="Helvetica", fontSize=8.9, leading=12.2, textColor=INK,
                         leftIndent=11, bulletIndent=3, spaceAfter=2.5),
 "ref":   ParagraphStyle("ref", fontName="Helvetica", fontSize=7.9, leading=10.6, textColor=INK,
                         leftIndent=11, firstLineIndent=-11, spaceAfter=2.5),
}

def inline(s):
    s = s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    s = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", s)
    s = re.sub(r"(?<!\*)\*([^*]+?)\*(?!\*)", r"<i>\1</i>", s)
    s = re.sub(r"`(.+?)`", r'<font face="Courier">\1</font>', s)
    s = s.replace("⟨", "&lt;").replace("⟩", "&gt;")
    return s

def md_table(block):
    rows = [r.strip().strip("|").split("|") for r in block if r.strip().startswith("|")]
    rows = [[c.strip() for c in r] for r in rows if not set("".join(r).replace(" ", "")) <= set("-:")]
    if not rows: return None
    ncol = max(len(r) for r in rows)
    rows = [r + [""]*(ncol-len(r)) for r in rows]
    cs = ParagraphStyle("c", fontName="Helvetica", fontSize=6.6, leading=8.4, textColor=INK)
    hs = ParagraphStyle("h", fontName="Helvetica-Bold", fontSize=6.6, leading=8.4, textColor=colors.white)
    data = [[Paragraph(inline(c), hs) for c in rows[0]]] + \
           [[Paragraph(inline(c), cs) for c in r] for r in rows[1:]]
    t = Table(data, colWidths=[BODY_W/ncol]*ncol, repeatRows=1, hAlign="LEFT")
    t.setStyle(TableStyle([
        ("BACKGROUND", (0,0), (-1,0), ACC),
        ("ROWBACKGROUNDS", (0,1), (-1,-1), [colors.white, colors.HexColor("#f6f6f2")]),
        ("GRID", (0,0), (-1,-1), 0.25, LINE),
        ("VALIGN", (0,0), (-1,-1), "MIDDLE"),
        ("LEFTPADDING", (0,0), (-1,-1), 3), ("RIGHTPADDING", (0,0), (-1,-1), 3),
        ("TOPPADDING", (0,0), (-1,-1), 2.2), ("BOTTOMPADDING", (0,0), (-1,-1), 2.2)]))
    return t

def fig_flowable(png, legend):
    iw, ih = PILImage.open(png).size
    w = BODY_W; h = w * ih / iw
    if h > 155*mm: h = 155*mm; w = h * iw / ih
    return KeepTogether([Image(png, width=w, height=h), Paragraph(inline(legend), S["cap"])])

FIGS = {
 1: ("Fig1_driver_detection.png",
     "**Figure 1 | Driver detection is strongest on ELN 2022 risk-defining lesions.** "
     "**a**, Per-driver AUROC under donor-grouped cross-validation, coloured by whether the lesion is "
     "risk-defining under ELN 2022. **b**, AUROC distribution for risk-defining versus other lesions "
     "(bars, means). **c**, Sensitivity versus specificity per driver; point area is proportional to the "
     "number of positive specimens; star marks the ideal corner."),
 2: ("Fig2_auroc_forest.png",
     "**Figure 2 | Every driver exceeds chance.** Per-driver AUROC with 95% donor-level bootstrap "
     "confidence intervals (B = 1,000); point area is proportional to positives, listed at right. "
     "All 28 drivers exceed a 1,000-shuffle label-permutation null (all P ≤ 0.004)."),
 3: ("Fig3_modality_and_augmentation.png",
     "**Figure 3 | The multimodal gain is not an artefact of RNA imputation, and bulk cohorts augment "
     "single-cell training.** **a**, Modality ablation ladder; removing all RNA-imputed blocks retains "
     "~70% of the AUROC gain over bulk RNA. **b**, Within BeatAML, imputed modalities improve a linear "
     "caller no more than a matched-width random nonlinear transform of the same RNA (n.s.). "
     "**c**, Pooling 707 bulk specimens into the training partition of the shared blocks improves the "
     "deployed model (Wilcoxon P = 0.016)."),
 4: ("Fig4_model_by_cohort.png",
     "**Figure 4 | Model and assay must be matched.** The bulk RNA-only caller performs well within bulk "
     "(BeatAML, Leucegene) but degrades markedly on single-cell input; the multimodal model, which "
     "requires single-cell modalities, is shown on its native assay. Dashed line separates assay types."),
 5: ("Fig5_calibration_utility_vaf.png",
     "**Figure 5 | Calibration, clinical utility and failure mode.** **a**, Reliability diagram on "
     "equal-count bins; nested Platt scaling reduces expected calibration error from 0.455 to 0.006. "
     "**b**, Decision-curve analysis: net benefit of the calibrated model against treat-all and "
     "treat-none references. **c**, Sensitivity stratified by variant allele frequency; missed lesions "
     "concentrate in subclonal disease."),
 6: ("Rx1_per_inhibitor_performance.png",
     "**Figure 6 | COMPASS-AML predicts ex-vivo inhibitor sensitivity to the limit of the assay.** "
     "Per-inhibitor performance across 118 agents in 520 specimens. The apparent AUROC of 0.774 falls "
     "to 0.671 when scored against the patient x drug interaction term alone, which is the only "
     "quantity constituting a drug-specific recommendation; 0.671 is ~92% of the directly measured "
     "assay reliability ceiling of 0.727. Predictability tracks assay reproducibility (Spearman 0.288, "
     "P = 0.0017) but not training-set size (0.112, P = 0.23)."),
 7: ("Sv1_survival_discrimination.png",
     "**Figure 7 | Survival discrimination by feature block against the clinical baseline.** Age and "
     "ELN risk are available at diagnosis and free, so the quantity of interest is the increment over "
     "them. The deployed model reaches C-index 0.756 +/- 0.029 across 60 re-draws of the sealed hold-out "
     "(the single sealed split gave 0.787) against 0.694 for age + "
     "ELN (+0.059, 95% CI +0.030 to +0.088). The molecular blocks alone do not beat the clinical "
     "baseline (-0.003); only the combination adds."),
 8: ("Sv4_tcga_external_validation.png",
     "**Figure 8 | External validation on TCGA-LAML with all coefficients frozen.** Kaplan-Meier by "
     "predicted risk tertile in 149 patients (92 deaths) profiled on a different platform. Cox "
     "coefficients, PCA rotation, gene selection and fusion weights were loaded unchanged; only the "
     "per-gene z-reference was cohort-matched. C-index 0.706 (95% CI 0.654-0.758); tertiles separate "
     "71.7% versus 13.7% two-year survival, log-rank P = 7.0 x 10^-10."),
}
TABLE_FILES = {
 1: ("Table1_pooled_heldout.tsv",
     "**Table 1 | Pooled held-out performance across three single-cell cohorts.** Every "
     "(specimen × driver) call with a known label is counted, with no minimum-positive filter. "
     "GSE281087 is scored panel-honestly: genes its panel never assayed are left unlabelled rather than "
     "counted as wild-type."),
}

def build():
    lines = open(SRC, encoding="utf-8").read().split("\n")
    story, i, placed = [], 0, set()
    while i < len(lines):
        L = lines[i]
        if L.startswith("| "):                                    # markdown table block
            blk = []
            while i < len(lines) and lines[i].strip().startswith("|"):
                blk.append(lines[i]); i += 1
            t = md_table(blk)
            if t: story += [Spacer(1, 3), t, Spacer(1, 7)]
            continue
        if L.startswith("# "):
            story.append(Paragraph(inline(L[2:]), S["title"]))
        elif L.startswith("### "):
            story.append(Paragraph(inline(L[4:]), S["h2"]))
        elif L.startswith("## "):
            head = L[3:].strip()
            if head in ("Display items",):                          # figures already inline
                i += 1
                while i < len(lines) and not lines[i].startswith("## "): i += 1
                continue
            story.append(Paragraph(inline(head), S["h1"]))
        elif L.startswith("> "):
            story.append(Paragraph(inline(L[2:]), S["quote"]))
        elif re.match(r"^\d+\. ", L):
            txt = L
            while i+1 < len(lines) and lines[i+1].startswith("   "):
                i += 1; txt += " " + lines[i].strip()
            story.append(Paragraph(inline(txt), S["ref"]))
        elif L.startswith("- "):
            story.append(Paragraph(inline(L[2:]), S["li"], bulletText="•"))
        elif L.strip() in ("---", ""):
            pass
        else:
            para = [L]
            while i+1 < len(lines) and lines[i+1].strip() and not re.match(
                    r"^(#|\||> |- |\d+\. |---)", lines[i+1]):
                i += 1; para.append(lines[i])
            txt = " ".join(para).strip()
            if txt:
                story.append(Paragraph(inline(txt), S["body"]))
                for n,(png,leg) in FIGS.items():                   # place a figure at its first callout
                    if n not in placed and re.search(r"\*\*Fig\. %d[a-c]?\*\*" % n, txt):
                        p = os.path.join(FIGD, png)
                        if os.path.exists(p):
                            story += [Spacer(1,4), fig_flowable(p, leg)]; placed.add(n)
                if 1 not in placed and "Table 1" in txt:
                    pass
        i += 1

    # any figure never cited inline goes at the end, plus Table 1
    for n,(png,leg) in FIGS.items():
        if n not in placed:
            p = os.path.join(FIGD, png)
            if os.path.exists(p): story += [Spacer(1,4), fig_flowable(p, leg)]
    story.append(PageBreak())
    story.append(Paragraph("Tables", S["h1"]))
    for n,(tsv,cap) in TABLE_FILES.items():
        f = os.path.join(TABD, tsv)
        if not os.path.exists(f): continue
        rows = list(csv.reader(open(f, encoding="utf-8"), delimiter="\t"))
        md = ["| " + " | ".join(rows[0]) + " |", "|" + "|".join(["-"]*len(rows[0])) + "|"] + \
             ["| " + " | ".join(r) + " |" for r in rows[1:] if any(r)]
        t = md_table(md)
        if t: story += [t, Paragraph(inline(cap), S["cap"])]
    for tsv in sorted(glob.glob(os.path.join(TABD, "SuppTable*.tsv"))):
        rows = list(csv.reader(open(tsv, encoding="utf-8"), delimiter="\t"))
        name = os.path.basename(tsv).replace(".tsv","").replace("_"," ")
        md = ["| " + " | ".join(rows[0]) + " |", "|" + "|".join(["-"]*len(rows[0])) + "|"] + \
             ["| " + " | ".join(r) + " |" for r in rows[1:] if any(r)]
        t = md_table(md)
        if t:
            story += [Spacer(1,6), Paragraph(inline("**%s**" % name), S["h2"]), t]

    def page(canv, doc):
        canv.saveState()
        canv.setFont("Helvetica", 7); canv.setFillColor(MUT)
        canv.drawString(20*mm, 12*mm, "MOSAIC-AML — manuscript draft")
        canv.drawRightString(A4[0]-20*mm, 12*mm, "page %d" % doc.page)
        canv.setStrokeColor(LINE); canv.setLineWidth(0.4)
        canv.line(20*mm, 15*mm, A4[0]-20*mm, 15*mm)
        canv.restoreState()

    doc = BaseDocTemplate(OUT, pagesize=A4, leftMargin=20*mm, rightMargin=20*mm,
                          topMargin=17*mm, bottomMargin=19*mm, title="MOSAIC-AML manuscript")
    doc.addPageTemplates([PageTemplate(id="n",
        frames=[Frame(20*mm, 19*mm, BODY_W, A4[1]-36*mm, id="f")], onPage=page)])
    doc.build(story)
    print("wrote %s (%.1f KB)" % (OUT, os.path.getsize(OUT)/1024))
    print("figures embedded:", sorted(placed) or "none inline")

build()
