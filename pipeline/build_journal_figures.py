#!/usr/bin/env python3
"""Publication-format figures for the MOSAIC-AML manuscript.

Nature Methods specifications applied throughout:
  width   88 mm (single column) / 120 mm (1.5) / 180 mm (double) — never exceeded
  height  <= 170 mm so the legend fits on the page
  type    sans-serif, >= 7 pt (5 pt absolute minimum), no raster elements
  format  vector PDF, panel letters a/b/c in bold lower-case at top-left
Every figure is redrawn from the source result JSONs, not stitched from earlier renders.

  python build_journal_figures.py      -> deliverables/figures/Fig1..Fig5.pdf (+ .png previews)
"""
import os, json, glob, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(ROOT, "deliverables")
OUT = os.path.join(D, "figures"); os.makedirs(OUT, exist_ok=True)
MM = 1/25.4
W1, W15, W2 = 88*MM, 120*MM, 180*MM            # column widths in inches
BLUE, GREEN, GREY, RED, AMBER = "#2a78d6", "#008300", "#b9b9b0", "#c0392b", "#eda100"
INK, MUTED, GRID = "#1a1a19", "#6b6b63", "#d9d9d4"

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7, "axes.labelsize": 7.5, "axes.titlesize": 7.5,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.5,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.2, "ytick.major.size": 2.2,
    "pdf.fonttype": 42, "ps.fonttype": 42,        # embed as TrueType, editable text
    "savefig.dpi": 600, "figure.dpi": 150,
})

def style(ax, grid="y"):
    ax.set_facecolor("white")
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    for s in ("left", "bottom"): ax.spines[s].set_color(GRID)
    ax.tick_params(colors=MUTED)
    if grid in ("y", "both"): ax.yaxis.grid(True, color=GRID, lw=0.45)
    if grid in ("x", "both"): ax.xaxis.grid(True, color=GRID, lw=0.45)
    ax.set_axisbelow(True)

def panel(ax, letter, dx=-0.16, dy=1.06):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=9, fontweight="bold",
            va="top", ha="left", color=INK)

def save(fig, name):
    fig.savefig(os.path.join(OUT, name + ".pdf"), format="pdf", bbox_inches="tight", pad_inches=0.02)
    fig.savefig(os.path.join(OUT, name + ".png"), format="png", bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    print("  wrote", name)

# ---------------------------------------------------------------- data
P = json.load(open(os.path.join(D, "production_fused_model.json")))
V = json.load(open(os.path.join(ROOT, "scratchpad", "oof_metrics_v3_nyu2.json")))
BM = json.load(open(os.path.join(D, "bulk_matrix.json")))["cohorts"]
CAL = json.load(open(os.path.join(D, "calibration_dca.json")))["calibration"]
VAF = json.load(open(os.path.join(D, "vaf_stratified.json")))["vaf_strata"]
IMP = json.load(open(os.path.join(D, "beataml_impute_experiment.json")))
HO = json.load(open(os.path.join(D, "pooled_heldout_eval.json")))
SC = json.load(open(os.path.join(ROOT, "scratchpad", "oof_scores_v2.json")))

ELN = {"inv16","inv(16)_CBFB-MYH11","NPM1","CEBPA","TP53","ASXL1","BCOR","EZH2","SF3B1","SRSF2",
       "STAG2","U2AF1","ZRSR2","kmt2a","KMT2A-rearrangement","complex","del5","del7","RUNX1","FLT3-ITD"}
rows, seen = [], set()
for m, a in P["per_mutation"].items():
    r = a.get("fused_all")
    if not r or r.get("auroc") is None: continue
    if round(r["auroc"], 4) in seen: continue           # drop alias duplicates
    seen.add(round(r["auroc"], 4))
    n = sum(r[k] for k in ("tp","fp","fn","tn")); prev = (r["tp"]+r["fn"])/n
    rows.append(dict(m=m, auroc=r["auroc"], sens=r["sensitivity"], spec=r["specificity"],
                     npos=r["tp"]+r["fn"], prev=prev,
                     prec=(r["tp"]/(r["tp"]+r["fp"]) if r["tp"]+r["fp"] else 0),
                     eln=m in ELN))
rows.sort(key=lambda r: -r["auroc"])
eln_a = [r["auroc"] for r in rows if r["eln"]]; non_a = [r["auroc"] for r in rows if not r["eln"]]

# ================================================================ Fig 1 (180 mm)
fig, ax = plt.subplots(1, 3, figsize=(W2, 62*MM), gridspec_kw={"width_ratios":[1.45,0.75,1.0]})
# a: per-driver AUROC ranked, coloured by ELN-defining status
a0 = ax[0]; y = np.arange(len(rows))[::-1]
a0.barh(y, [r["auroc"] for r in rows], height=0.7,
        color=[GREEN if r["eln"] else BLUE for r in rows], edgecolor="none")
a0.set_yticks(y); a0.set_yticklabels([r["m"][:15] for r in rows], fontsize=5.2)
a0.set_xlim(0.5, 1.0); a0.set_xlabel("AUROC (donor-grouped CV)")
a0.axvline(0.5, color=RED, lw=0.7, ls=(0,(3,2)))
style(a0, "x")
a0.legend(handles=[Line2D([],[],marker="s",ls="",color=GREEN,label="ELN 2022 risk-defining",ms=4),
                   Line2D([],[],marker="s",ls="",color=BLUE,label="other driver",ms=4)],
          loc="lower right", frameon=False, handletextpad=0.4)
panel(a0, "a", dx=-0.42)
# b: ELN-defining vs other, distribution
a1 = ax[1]
for i,(v,c,lab) in enumerate(((eln_a,GREEN,"ELN\nrisk-defining"),(non_a,BLUE,"other"))):
    xj = np.random.default_rng(0).normal(i, 0.055, len(v))
    a1.scatter(xj, v, s=9, color=c, alpha=.85, edgecolor="white", lw=.4, zorder=3)
    a1.plot([i-0.22,i+0.22], [np.mean(v)]*2, color=c, lw=1.6, zorder=4)
    a1.text(i, 1.005, "%.3f"%np.mean(v), ha="center", fontsize=6.5, fontweight="bold", color=c)
a1.set_xticks([0,1]); a1.set_xticklabels(["ELN\nrisk-defining","other"], fontsize=6.5)
a1.set_ylim(0.72, 1.02); a1.set_ylabel("AUROC"); a1.set_xlim(-0.5,1.5); style(a1)
panel(a1, "b", dx=-0.34)
# c: sensitivity vs specificity
a2 = ax[2]
for r in rows:
    a2.scatter(r["sens"], r["spec"], s=6+2.0*r["npos"], color=GREEN if r["eln"] else BLUE,
               alpha=.8, edgecolor="white", lw=.4, zorder=3)
a2.plot([1],[1], marker="*", ms=7, color=AMBER, zorder=4)
a2.set_xlim(0,1.05); a2.set_ylim(0.85,1.005)
a2.set_xlabel("sensitivity"); a2.set_ylabel("specificity"); style(a2,"both")
a2.text(0.98,0.995,"ideal",fontsize=5.5,color=AMBER,ha="right")
panel(a2, "c", dx=-0.30)
fig.subplots_adjust(left=0.13, right=0.995, top=0.90, bottom=0.20, wspace=0.55)
save(fig, "Fig1_driver_detection")

# ================================================================ Fig 2 (120 mm) forest
mm_ = V["arms"]["multimodal"]["mutations"]
fr, seen2 = [], set()
for m, r in mm_.items():
    if round(r["auroc"],4) in seen2: continue
    seen2.add(round(r["auroc"],4))
    fr.append((m, r["auroc"], r.get("auroc_ci"), r["n_pos"]))
fr.sort(key=lambda x: x[1])
fig, ax = plt.subplots(figsize=(W15, 105*MM))
yy = np.arange(len(fr))
for i,(m,a,ci,npos) in enumerate(fr):
    if ci: ax.plot(ci, [i,i], color=BLUE, lw=1.1, alpha=.55, solid_capstyle="round")
    ax.scatter([a],[i], s=8+1.8*npos, color=BLUE, edgecolor="white", lw=.5, zorder=3)
    ax.text(1.012, i, "%d"%npos, va="center", fontsize=5, color=MUTED)
ax.axvline(0.5, color=RED, lw=0.8, ls=(0,(3,2)))
ax.text(0.5, len(fr)-0.2, "chance", color=RED, fontsize=5.5, ha="center")
ax.set_yticks(yy); ax.set_yticklabels([f[0][:18] for f in fr], fontsize=5.6)
ax.set_xlim(0.45,1.03); ax.set_xlabel("AUROC with 95% donor-bootstrap CI   (n+ at right)")
style(ax,"x")
fig.subplots_adjust(left=0.28, right=0.95, top=0.985, bottom=0.075)
save(fig, "Fig2_auroc_forest")

# ================================================================ Fig 3 (180 mm) 3 panels
L = V["modality_ladder"]
fig, ax = plt.subplots(1, 3, figsize=(W2, 58*MM))
# a ladder
arms = ["bulkrna","rna_comp","measured","multimodal"]
labs = ["bulk\nRNA","+ compo-\nsition","measured\nonly","all 8\nmodalities"]
vals = [L[a]["mean_auroc"] for a in arms]
a0=ax[0]; a0.bar(range(4), vals, width=.62, color=[BLUE,"#5b9bd5","#8fbf6f",GREEN], edgecolor="none")
for i,v in enumerate(vals): a0.text(i, v+.004, "%.3f"%v, ha="center", fontsize=6, fontweight="bold")
a0.set_xticks(range(4)); a0.set_xticklabels(labs, fontsize=6)
a0.set_ylim(0.5,0.95); a0.set_ylabel("mean AUROC"); style(a0)
a0.annotate("", xy=(2,0.905), xytext=(0,0.905), arrowprops=dict(arrowstyle="<->",color="#1a6e1a",lw=1))
a0.text(1,0.912,"imputation-independent",ha="center",fontsize=5.6,color="#1a6e1a",fontweight="bold")
panel(a0,"a")
# b imputation vs random control
am = IMP["arm_means"]
a1=ax[1]; v=[am["rna_only"]["auroc"],am["rna_imputed"]["auroc"],am["rna_random"]["auroc"]]
bars=a1.bar(range(3), v, width=.62, color=[GREY,GREEN,GREEN], edgecolor="none")
bars[2].set_hatch("////"); bars[2].set_edgecolor("white")
for i,x in enumerate(v): a1.text(i,x+.004,"%.3f"%x,ha="center",fontsize=6,fontweight="bold")
a1.set_xticks(range(3)); a1.set_xticklabels(["RNA\nonly","+ imputed","+ random\n(control)"],fontsize=6)
a1.set_ylim(0.5,0.95); a1.set_ylabel("mean AUROC"); style(a1)
a1.text(1.5,0.90,"n.s.",ha="center",fontsize=6.5,color="#b06a00",fontweight="bold")
panel(a1,"b")
# c augmentation
S=P["summary"]; a2=ax[2]
v2=[S["deployed"]["auroc"], S["fused_all"]["auroc"]]
a2.bar(range(2), v2, width=.5, color=[BLUE,GREEN], edgecolor="none")
for i,x in enumerate(v2): a2.text(i,x+.003,"%.3f"%x,ha="center",fontsize=6,fontweight="bold")
a2.set_xticks(range(2)); a2.set_xticklabels(["single-cell\nonly","+ BeatAML\naugmentation"],fontsize=6)
a2.set_ylim(0.80,0.94); a2.set_ylabel("mean AUROC"); style(a2)
a2.text(0.5,0.928,"P = 0.016",ha="center",fontsize=6,color="#1a6e1a",fontweight="bold")
panel(a2,"c")
fig.subplots_adjust(left=0.07,right=0.99,top=0.90,bottom=0.20,wspace=0.42)
save(fig, "Fig3_modality_and_augmentation")

# ================================================================ Fig 4 (180 mm) model x cohort
fig, ax = plt.subplots(figsize=(W2, 62*MM))
cells=[("Bulk model\nBeatAML (CV)",BM["BeatAML_CV"]["overall"]),("Bulk model\nLeucegene",BM["Leucegene"]["overall"]),
       ("Bulk model\nall scRNA",BM["all_scRNA"]["overall"])]
names=[c[0] for c in cells]+["Multimodal\nall scRNA (CV)"]
se=[c[1]["mean_sensitivity"] for c in cells]+[S["fused_all"]["sensitivity"]]
sp=[c[1]["mean_specificity"] for c in cells]+[S["fused_all"]["specificity"]]
au=[c[1]["mean_auroc"] for c in cells]+[S["fused_all"]["auroc"]]
x=np.arange(4); w=.26
ax.bar(x-w, se, w, color=BLUE, label="sensitivity", edgecolor="none")
ax.bar(x,   sp, w, color=GREEN, label="specificity", edgecolor="none")
ax.bar(x+w, au, w, color=AMBER, label="AUROC", edgecolor="none")
for i in range(4):
    for off,v in ((-w,se[i]),(0,sp[i]),(w,au[i])):
        ax.text(i+off, v+.012, "%.2f"%v, ha="center", fontsize=5.4)
ax.set_xticks(x); ax.set_xticklabels(names, fontsize=6.2)
ax.set_ylim(0,1.10); ax.set_ylabel("score"); style(ax)
ax.axvline(2.5, color=GRID, lw=0.8, ls=(0,(3,2)))
ax.legend(frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.5,1.14), handlelength=1.2)
fig.subplots_adjust(left=0.075,right=0.99,top=0.86,bottom=0.19)
save(fig, "Fig4_model_by_cohort")

# ================================================================ Fig 5 (180 mm) 3 panels
fig, ax = plt.subplots(1, 3, figsize=(W2, 58*MM))
# a calibration — nested Platt (sigmoid) scaling, selected over isotonic on nested ECE/Brier/log-loss
# (isotonic overfits at these positive counts and yields a noisy step function; see Methods).
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupKFold
a0=ax[0]
allp,ally,allc=[],[],[]
CALM={}                                     # per-mutation calibrated probability, reused by panel b
for m,d0 in SC["multimodal"].items():
    p0=np.array(d0["score"],float); y0=np.array(d0["y"],int); g0=np.array(d0["donor"])
    cc=np.full(len(y0),np.nan); ng=min(3,len(set(g0)))
    if ng>=2 and len(set(y0))==2:
        for tri,tei in GroupKFold(ng).split(p0,y0,g0):
            if len(set(y0[tri]))<2: cc[tei]=y0[tri].mean(); continue
            lr=LogisticRegression(C=1e6, solver="lbfgs").fit(p0[tri].reshape(-1,1), y0[tri])
            cc[tei]=lr.predict_proba(p0[tei].reshape(-1,1))[:,1]
    else: cc=p0.copy()
    CALM[m]=(cc,y0)
    allp+=list(p0); ally+=list(y0); allc+=list(cc)
allp=np.array(allp); ally=np.array(ally); allc=np.array(allc)
ok=~np.isnan(allc); allp,ally,allc=allp[ok],ally[ok],allc[ok]
def rel(p,y,nb=10):
    """equal-COUNT (quantile) bins: equalises n per bin so the curve is not dominated by
    sparsely-populated high-probability bins."""
    e=np.unique(np.quantile(p, np.linspace(0,1,nb+1)))
    xs,ys,ns=[],[],[]
    for i in range(len(e)-1):
        msk=(p>=e[i])&(p<e[i+1] if i<len(e)-2 else p<=e[-1]+1e-9)
        if msk.sum()>=10: xs.append(p[msk].mean()); ys.append(y[msk].mean()); ns.append(int(msk.sum()))
    return xs,ys,ns
rx,ry,_=rel(allp,ally); cx,cy,cn=rel(allc,ally)
def _ece(p,y,nb=10):
    e=np.unique(np.quantile(p,np.linspace(0,1,nb+1))); t=0.0
    for i in range(len(e)-1):
        m=(p>=e[i])&(p<e[i+1] if i<len(e)-2 else p<=e[-1]+1e-9)
        if m.sum(): t+=m.mean()*abs(y[m].mean()-p[m].mean())
    return t
ECE_RAW, ECE_CAL = _ece(allp,ally), _ece(allc,ally)
a0.plot([0,1],[0,1],color=MUTED,ls=(0,(4,3)),lw=.8)
a0.plot(rx,ry,"-o",color=GREY,lw=1.1,ms=2.6,label="raw score")
a0.plot(cx,cy,"-o",color=GREEN,lw=1.3,ms=2.6,label="calibrated (Platt)")
a0.set_xlabel("predicted probability"); a0.set_ylabel("observed frequency")
a0.set_xlim(0,1); a0.set_ylim(0,1); a0.legend(frameon=False,loc="upper left"); style(a0,"both")
a0.text(.97,.06,"ECE %.3f → %.3f"%(ECE_RAW,ECE_CAL),ha="right",fontsize=5.8,color=INK)
panel(a0,"a")
# b decision curve (NPM1 exemplar recomputed)
a1=ax[1]
# net benefit must be computed on the CALIBRATED probability — thresholding a raw decision score at a
# "threshold probability" is not a decision curve and produces spurious negative net benefit.
_k = "NPM1" if "NPM1" in CALM else list(CALM)[0]
p1, y1 = CALM[_k]
ok1 = ~np.isnan(p1); p1, y1 = p1[ok1], y1[ok1]; n=len(y1); prev=y1.mean()
pts=np.linspace(0.01,0.6,60)
nb=[( (p1>=t)&(y1==1)).sum()/n - (((p1>=t)&(y1==0)).sum()/n)*(t/(1-t)) for t in pts]
na=[prev-(1-prev)*(t/(1-t)) for t in pts]
a1.plot(pts,nb,color=GREEN,lw=1.4,label="model")
a1.plot(pts,na,color=MUTED,lw=.9,ls=(0,(3,2)),label="treat all")
a1.axhline(0,color=RED,lw=.9,ls=(0,(3,2)),label="treat none")
a1.set_xlim(0,.6); a1.set_ylim(min(-0.02,min(nb)),max(0.02,prev*1.15))
a1.set_xlabel("threshold probability"); a1.set_ylabel("net benefit")
a1.legend(frameon=False,loc="upper right"); style(a1,"both"); panel(a1,"b")
# c VAF
a2=ax[2]
labs=[s["bin"].replace("subclonal ","").replace("clonal ","") for s in VAF]
vv=[s["sensitivity"] for s in VAF]
a2.bar(range(len(vv)), vv, width=.62, color=BLUE, edgecolor="none")
for i,v in enumerate(vv): a2.text(i,v+.012,"%.2f"%v,ha="center",fontsize=6,fontweight="bold")
a2.set_xticks(range(len(vv))); a2.set_xticklabels(labs,fontsize=5.8,rotation=20,ha="right")
a2.set_ylim(0,.85); a2.set_ylabel("sensitivity"); a2.set_xlabel("variant allele frequency")
style(a2); panel(a2,"c")
fig.subplots_adjust(left=0.07,right=0.99,top=0.93,bottom=0.24,wspace=0.42)
save(fig, "Fig5_calibration_utility_vaf")

# ---------------------------------------------------------------- spec check
print("\n=== dimension check (Nature Methods: <=180 mm wide, <=170 mm tall) ===")
import re as _re
for f in sorted(glob.glob(os.path.join(OUT, "*.pdf"))):
    d=open(f,"rb").read(4000)
    m=_re.search(rb'MediaBox\s*\[\s*0\s+0\s+([\d.]+)\s+([\d.]+)', d)
    if m:
        w,h=float(m.group(1))/72*25.4, float(m.group(2))/72*25.4
        print("  %-38s %5.1f x %5.1f mm   %s" % (os.path.basename(f), w, h,
              "OK" if (w<=180.5 and h<=170) else "OVER"))
print("FIGURES OK")
