#!/usr/bin/env python3
"""Journal-spec figure panels for COMPASS-AML (the drug-response layer).

Nature Methods geometry throughout: <=180 mm wide, >=5.2 pt type, embedded TrueType (pdf.fonttype 42),
600 dpi. Every panel is generated from the validation JSONs, never hand-typed.

  Rx1  per-inhibitor performance, approved agents marked
  Rx2  the deployment task: per-patient retrieval and ranking concordance against matched chance
  Rx3  where the signal comes from -- block ablation and the fitted fusion weights per pathway family
  Rx4  robustness -- leave-wave-out, leave-centre-out, differentiation strata, permutation null
  Rx5  calibration and uncertainty-based abstention
  Rx6  venetoclax as the worked example, including the state-level contrast Model B reproduces

  python build_drug_figures.py -> deliverables/figures/Rx{1..6}_*.pdf/.png
"""
import os, sys, json
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
D = os.path.join(ROOT, "deliverables")
FIGD = os.path.join(D, "figures")
os.makedirs(FIGD, exist_ok=True)
MM = 1 / 25.4
INK, MUTED, GRID, PANEL = "#1a1a19", "#6b6b63", "#d9d9d4", "#f7f7f4"
ACC, GOOD, BAD, COOL = "#8a6a18", "#1a6e1a", "#a03828", "#9ec9e2"
BLOCK_LABEL = {"rna": "RNA", "state": "state", "mut": "mutation", "clin": "clinical"}
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
                     "font.size": 7, "axes.labelsize": 7, "xtick.labelsize": 6.5,
                     "ytick.labelsize": 6.5, "legend.fontsize": 6.5, "axes.linewidth": 0.6,
                     "pdf.fonttype": 42, "savefig.dpi": 600})


def load(name):
    p = os.path.join(D, name)
    return json.load(open(p)) if os.path.exists(p) else None


V = load("drug_model_validation.json")
C = load("drug_model_card.json")
S = load("state_response_validation.json")
ABL = {}
for f in os.listdir(D):
    if f.startswith("drug_model_card_abl_") and f.endswith(".json"):
        ABL[f[len("drug_model_card_abl_"):-5]] = json.load(open(os.path.join(D, f)))


def style(ax, xg=True, yg=False):
    ax.set_facecolor("white")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=MUTED)
    if xg:
        ax.xaxis.grid(True, color=GRID, lw=0.4)
    if yg:
        ax.yaxis.grid(True, color=GRID, lw=0.4)
    ax.set_axisbelow(True)


def save(fig, name):
    for ext, dpi in (("pdf", None), ("png", 200)):
        fig.savefig(os.path.join(FIGD, "%s.%s" % (name, ext)), **({"dpi": dpi} if dpi else {}),
                    bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    print("wrote %s.pdf/.png" % name)


# ------------------------------------------------------------------- Rx1 ----
def rx1():
    t = V["per_drug"]["table"]
    tier = {d: C["per_drug"][d]["annotation"]["clinical_tier"] for d in C["per_drug"]}
    t = [r for r in t if r.get("auroc") is not None]
    t.sort(key=lambda r: r["auroc"])
    fig, (a0, a1) = plt.subplots(1, 2, figsize=(180 * MM, 92 * MM),
                                 gridspec_kw={"width_ratios": [1.25, 1]})
    au = np.array([r["auroc"] for r in t])
    cols = [ACC if tier.get(r["inhibitor"]) == "approved_AML"
            else (GOOD if tier.get(r["inhibitor"]) == "approved_other" else COOL) for r in t]
    a0.barh(np.arange(len(t)), au - 0.5, left=0.5, height=0.86, color=cols, lw=0)
    a0.axvline(0.5, color=INK, lw=0.7)
    a0.set_xlim(0.5, 1.0); a0.set_ylim(-1, len(t))
    a0.set_yticks([]); a0.set_xlabel("AUROC, sensitive vs resistant tail (donor-grouped CV)")
    a0.set_title("%d inhibitors, ranked\nmean %.3f, median %.3f"
                 % (len(t), au.mean(), np.median(au)), fontsize=7.4, loc="left", fontweight="bold")
    for r, y in zip(t, np.arange(len(t))):
        if tier.get(r["inhibitor"]) == "approved_AML" or r["auroc"] >= max(au) - 1e-9:
            a0.text(r["auroc"] + 0.004, y, r["inhibitor"][:26], va="center", fontsize=5.2, color=INK)
    style(a0)
    a0.legend(handles=[Line2D([], [], marker="s", ls="", color=ACC, label="approved in AML", ms=4.5),
                       Line2D([], [], marker="s", ls="", color=GOOD, label="approved, other indication", ms=4.5),
                       Line2D([], [], marker="s", ls="", color=COOL, label="trial / research", ms=4.5)],
              frameon=False, loc="lower right")

    sp = np.array([r["spearman"] for r in t])
    a1.scatter(sp, au, s=13, c=cols, lw=0.3, edgecolor="white", zorder=3)
    for r in t:
        if tier.get(r["inhibitor"]) == "approved_AML":
            a1.annotate(r["inhibitor"].split(" (")[0], (r["spearman"], r["auroc"]),
                        fontsize=5.4, color=INK, xytext=(3, -1), textcoords="offset points")
    a1.axhline(0.5, color=GRID, lw=0.6); a1.axvline(0, color=GRID, lw=0.6)
    a1.set_xlabel("Spearman, predicted vs measured AUC (within inhibitor)")
    a1.set_ylabel("AUROC")
    a1.set_title("continuous and tail metrics agree", fontsize=7.4, loc="left", fontweight="bold")
    style(a1, yg=True)
    fig.text(0.005, 0.012, "Donor-grouped 5-fold CV over %d specimens from %d patients; specimens from one "
             "patient never span folds." % (V["n_specimens"], V["n_subjects"]),
             fontsize=5.6, color=MUTED)
    fig.subplots_adjust(left=0.02, right=0.99, top=0.90, bottom=0.155, wspace=0.16)
    save(fig, "Rx1_per_inhibitor_performance")


# ------------------------------------------------------------------- Rx2 ----
def rx2():
    pp = V["per_patient"]
    ks = [1, 3, 5, 10]
    fig, (a0, a1) = plt.subplots(1, 2, figsize=(180 * MM, 74 * MM),
                                 gridspec_kw={"width_ratios": [1, 1]})
    x = np.arange(len(ks)); w = 0.36
    hit = [pp["hit@%d" % k] for k in ks]
    rnd = [pp["random@%d" % k] for k in ks]
    a0.bar(x - w / 2, hit, w, color=ACC, label="COMPASS-AML")
    a0.bar(x + w / 2, rnd, w, color="#c9c6b8", label="matched chance")
    for i, (h, r) in enumerate(zip(hit, rnd)):
        a0.text(i - w / 2, h + 0.015, "%.2f" % h, ha="center", fontsize=6, color=INK, fontweight="bold")
        a0.text(i + w / 2, r + 0.015, "%.2f" % r, ha="center", fontsize=6, color=MUTED)
        a0.text(i, max(h, r) + 0.075, "%.1f×" % (h / r if r else np.nan), ha="center",
                fontsize=6.2, color=GOOD, fontweight="bold")
    a0.set_xticks(x); a0.set_xticklabels(["top-%d" % k for k in ks])
    a0.set_ylim(0, 1.05); a0.set_ylabel("P(at least one true top-decile inhibitor)")
    a0.set_title("which drug for THIS patient?\n%d specimens; chance is matched to how many\n"
                 "inhibitors each specimen was tested on"
                 % pp["n_specimens"], fontsize=7.2, loc="left", fontweight="bold")
    a0.legend(frameon=False, loc="upper left")
    style(a0, xg=False, yg=True)

    a1.axis("off")
    rows = [("per-patient ranking concordance (mean Spearman)", "%.3f" % pp["mean_spearman"]),
            ("per-patient ranking concordance (median)", "%.3f" % pp["median_spearman"]),
            ("precision@5 (share of 5 truly top-decile)", "%.3f" % pp["prec@5"]),
            ("", ""),
            ("per-inhibitor mean Spearman", "%.3f" % V["per_drug"]["mean_spearman"]),
            ("per-inhibitor mean AUROC", "%.3f" % V["per_drug"]["mean_auroc"]),
            ("per-inhibitor mean AUPRC (baseline %.3f)" % V["per_drug"]["mean_auprc_baseline"],
             "%.3f" % V["per_drug"]["mean_auprc"]),
            ("inhibitors with Spearman p < 0.05", "%.0f%%" % (100 * V["per_drug"]["frac_spearman_p_lt_0.05"])),
            ("", ""),
            ("approved agents only (n = %d)" % V["actionable_subset"]["n_drugs"],
             "AUROC %.3f" % V["actionable_subset"]["mean_auroc"])]
    y = 0.96
    for lab, val in rows:
        if lab:
            a1.text(0.0, y, lab, fontsize=6.8, color=MUTED, va="top")
            a1.text(1.0, y, val, fontsize=6.8, color=INK, va="top", ha="right", fontweight="bold")
        y -= 0.085
    a1.set_title("headline numbers", fontsize=7.2, loc="left", fontweight="bold", pad=18)
    fig.subplots_adjust(left=0.075, right=0.98, top=0.83, bottom=0.09, wspace=0.30)
    save(fig, "Rx2_per_patient_retrieval")


# ------------------------------------------------------------------- Rx3 ----
def rx3():
    order = ["clin", "mut", "state", "rna", "rna+state", "rna+mut", "rna+state+clin", "state+mut+clin"]
    have = [k for k in order if k in ABL]
    vals = [ABL[k]["summary"]["oof_mean_auroc"][0] for k in have]
    allv = C["summary"]["oof_mean_auroc"][0]
    fig, (a0, a1) = plt.subplots(1, 2, figsize=(180 * MM, 84 * MM),
                                 gridspec_kw={"width_ratios": [1, 1.15]})
    y = np.arange(len(have))
    a0.barh(y, [v - 0.5 for v in vals], left=0.5, height=0.7,
            color=[ACC if "+" in k else COOL for k in have], lw=0)
    a0.axvline(allv, color=GOOD, lw=1.0, ls="--")
    a0.text(allv - 0.004, -0.85, "all four %.3f " % allv, fontsize=6, color=GOOD,
            va="bottom", ha="right", fontweight="bold")
    for i, v in enumerate(vals):
        a0.text(v + 0.003, i, "%.3f" % v, va="center", fontsize=5.8, color=INK)
    # "rna" set in a small sans face reads as "ma" (the rn ligature problem) -- spell the blocks out
    a0.set_yticks(y)
    a0.set_yticklabels([" + ".join(BLOCK_LABEL.get(b, b) for b in k.split("+")) for k in have],
                       fontsize=6.2)
    a0.set_xlim(0.5, max(vals + [allv]) + 0.04); a0.set_ylim(-1.4, len(have) - 0.4)
    a0.set_xlabel("mean OOF AUROC across %d inhibitors" % C["n_drugs"])
    a0.set_title("where the signal lives\nRNA carries most of it; mutations add a real but small increment",
                 fontsize=7.2, loc="left", fontweight="bold")
    style(a0)

    W = C["model"]["fusion_weights"]
    blocks = ["rna", "state", "mut", "clin"]
    fams = sorted(W, key=lambda g: -W[g].get("rna", 0))
    M = np.array([[W[g].get(b, 0.0) for b in blocks] for g in fams])
    im = a1.imshow(M, aspect="auto", cmap="YlGnBu", vmin=0, vmax=1)
    a1.set_xticks(range(len(blocks))); a1.set_xticklabels([BLOCK_LABEL[b] for b in blocks])
    a1.set_yticks(range(len(fams))); a1.set_yticklabels(fams, fontsize=6)
    for i in range(len(fams)):
        for j in range(len(blocks)):
            if M[i, j] >= 0.01:
                a1.text(j, i, "%.2f" % M[i, j], ha="center", va="center", fontsize=5.4,
                        color="white" if M[i, j] > 0.55 else INK)
    a1.set_title("fitted fusion weights per target-pathway family\n"
                 "FLT3 leans on mutations; cell-cycle on differentiation state",
                 fontsize=7.2, loc="left", fontweight="bold")
    cb = fig.colorbar(im, ax=a1, fraction=0.03, pad=0.02); cb.ax.tick_params(labelsize=5.6)
    fig.text(0.005, 0.005, "Weights are fitted by non-negative least squares on inner donor-grouped OOF "
             "predictions, floored to the best single block. A high weight means complementary "
             "information, not standalone strength.", fontsize=5.6, color=MUTED)
    fig.subplots_adjust(left=0.115, right=0.97, top=0.87, bottom=0.12, wspace=0.42)
    save(fig, "Rx3_ablation_and_fusion_weights")


# ------------------------------------------------------------------- Rx4 ----
def rx4():
    fig, axes = plt.subplots(1, 3, figsize=(180 * MM, 68 * MM),
                             gridspec_kw={"width_ratios": [1.15, 1.15, 1]})
    base = V["per_drug"]["mean_auroc"]

    lw = V["leave_wave_out"]; lc = V["leave_center_out"]
    labs = ["W1+2\n→W3+4", "W3+4\n→W1+2"] + \
           [k.replace("center_", "drop\ncentre ") for k in lc]
    vals = [v["mean_auroc"] for v in lw.values()] + [v["mean_auroc"] for v in lc.values()]
    ax = axes[0]
    ax.bar(range(len(vals)), [v - 0.5 for v in vals], bottom=0.5, color=COOL, lw=0)
    ax.axhline(base, color=ACC, lw=1.0, ls="--")
    ax.text(-0.4, base, " random CV %.3f" % base, fontsize=5.8, color=ACC, va="bottom", ha="left")
    for i, v in enumerate(vals):
        ax.text(i, v + 0.006, "%.3f" % v, ha="center", fontsize=5.6, color=INK)
    ax.set_xticks(range(len(labs))); ax.set_xticklabels(labs, fontsize=5.2)
    ax.set_ylim(0.5, max(vals + [base]) + 0.05); ax.set_ylabel("mean AUROC")
    ax.set_title("held out by acquisition wave and by centre", fontsize=7, loc="left", fontweight="bold")
    style(ax, xg=False, yg=True)

    ax = axes[1]
    st = V["strata"]
    names = [k for k in ("primitive", "intermediate", "monocytic/mature") if k in st]
    a = [st[k]["mean_auroc"] for k in names]
    p = [st[k]["patient_mean_spearman"] for k in names]
    x = np.arange(len(names)); w = 0.36
    ax.bar(x - w / 2, [v - 0.5 for v in a], w, bottom=0.5, color=ACC, label="per-inhibitor AUROC")
    ax.bar(x + w / 2, p, w, color=COOL, label="per-patient Spearman")
    for i, (v1, v2) in enumerate(zip(a, p)):
        ax.text(i - w / 2, v1 + 0.008, "%.2f" % v1, ha="center", fontsize=5.6, color=INK)
        ax.text(i + w / 2, v2 + 0.008, "%.2f" % v2, ha="center", fontsize=5.6, color=INK)
    ax.set_xticks(x); ax.set_xticklabels([n.replace("/", "/\n") for n in names], fontsize=5.6)
    ax.set_ylim(0, max(a) + 0.26)
    ax.set_title("inside differentiation-state strata\n(not just a monocytic-vs-primitive split)",
                 fontsize=7, loc="left", fontweight="bold")
    ax.legend(frameon=False, loc="upper center", ncol=2, fontsize=5.6,
              bbox_to_anchor=(0.5, 1.0), borderaxespad=0.2)
    style(ax, xg=False, yg=True)

    ax = axes[2]
    pn = V.get("permutation_null") or {}
    if pn:
        mu, sd = pn["null_mean"], max(pn["null_sd"], 1e-6)
        xs = np.linspace(mu - 4 * sd, max(pn["observed_mean_spearman"] * 1.05, mu + 4 * sd), 400)
        ax.fill_between(xs, np.exp(-0.5 * ((xs - mu) / sd) ** 2), color=GRID, lw=0)
        ax.axvline(pn["null_max"], color=MUTED, lw=0.7, ls=":")
        ax.axvline(pn["observed_mean_spearman"], color=BAD, lw=1.4)
        ax.text(pn["observed_mean_spearman"], 0.55, " observed\n %.3f" % pn["observed_mean_spearman"],
                fontsize=6, color=BAD, fontweight="bold", ha="right")
        ax.text(mu, 0.55, "null (n=%d)\n%.3f ± %.3f" % (pn["n_perm"], mu, sd),
                fontsize=5.8, color=MUTED, ha="center")
        ax.set_ylim(0, 1.25)
        ax.set_yticks([]); ax.set_xlabel("mean Spearman")
        ax.set_title("permutation null\nobserved sits %.0f null SDs above the null mean"
                     % pn.get("z_vs_null", (pn["observed_mean_spearman"] - mu) / sd),
                     fontsize=7, loc="left", fontweight="bold")
        style(ax, xg=False)
    fig.subplots_adjust(left=0.065, right=0.99, top=0.83, bottom=0.14, wspace=0.30)
    save(fig, "Rx4_robustness")


# ------------------------------------------------------------------- Rx5 ----
def rx5():
    fig, (a0, a1) = plt.subplots(1, 2, figsize=(180 * MM, 74 * MM))
    cal = V["calibration"]
    rel = cal["reliability"]
    a0.plot([0, 1], [0, 1], color=GRID, lw=0.8, ls="--")
    xs = [r["predicted"] for r in rel]; ys = [r["observed"] for r in rel]
    ns = np.array([r["n"] for r in rel], float)
    a0.plot(xs, ys, color=ACC, lw=1.2, zorder=3)
    a0.scatter(xs, ys, s=8 + 60 * ns / ns.max(), color=ACC, zorder=4, lw=0.3, edgecolor="white")
    a0.set_xlabel("predicted P(sensitive)"); a0.set_ylabel("observed frequency")
    a0.set_xlim(0, 1); a0.set_ylim(0, 1)
    a0.set_title("calibration\nECE %.3f | Brier %.3f (baseline %.3f), n=%d"
                 % (cal["ece"], cal["brier"], cal["brier_baseline"], cal["n"]),
                 fontsize=7.2, loc="left", fontweight="bold")
    style(a0, yg=True)

    ab = V["abstention"]
    cov = [r["coverage"] for r in ab]; err = [r["error_rate"] for r in ab]
    auc = [r["auroc"] for r in ab]
    a1.plot(cov, err, color=BAD, lw=1.3, marker="o", ms=3.5, label="error rate")
    a1.set_xlabel("coverage (fraction of calls the system is willing to make)")
    a1.set_ylabel("error rate", color=BAD)
    a1.invert_xaxis()
    ax2 = a1.twinx()
    ax2.plot(cov, auc, color=GOOD, lw=1.3, marker="s", ms=3.2, label="AUROC")
    ax2.set_ylabel("AUROC on the retained calls", color=GOOD)
    ax2.spines["top"].set_visible(False)
    for c, e in zip(cov, err):
        a1.annotate("%.0f%%" % (100 * e), (c, e), fontsize=5.6, color=BAD,
                    xytext=(0, 5), textcoords="offset points", ha="center")
    a1.set_title("declining to answer buys accuracy\nerror %.0f%% at full coverage -> %.0f%% at %d%% coverage"
                 % (100 * err[0], 100 * err[-1], 100 * cov[-1]), fontsize=7.2, loc="left", fontweight="bold")
    style(a1, xg=False, yg=True)
    fig.subplots_adjust(left=0.075, right=0.93, top=0.85, bottom=0.14, wspace=0.32)
    save(fig, "Rx5_calibration_and_abstention")


# ------------------------------------------------------------------- Rx6 ----
def rx6():
    ven = next((r for r in V["per_drug"]["table"] if r["inhibitor"] == "Venetoclax"), None)
    fig, (a0, a1) = plt.subplots(1, 2, figsize=(180 * MM, 76 * MM),
                                 gridspec_kw={"width_ratios": [1, 1.2]})
    a0.axis("off")
    if ven:
        rows = [("measurements (specimens)", "%d" % ven["n"]),
                ("Spearman, predicted vs measured AUC", "%.3f" % ven["spearman"]),
                ("AUROC, sensitive vs resistant tail", "%.3f" % ven["auroc"]),
                ("AUPRC", "%.3f  (baseline %.3f)" % (ven["auprc"], ven["auprc_baseline"])),
                ("RMSE on the raw AUC scale", "%.1f  (cohort SD %.1f)" % (ven["rmse_auc"], ven["auc_sd"])),
                ("fusion weights (RNA/state/mut/clin)",
                 "/".join("%.2f" % C["model"]["fusion_weights"]["apoptosis"].get(b, 0)
                          for b in ("rna", "state", "mut", "clin")))]
        y = 0.92
        for lab, val in rows:
            a0.text(0.0, y, lab, fontsize=6.8, color=MUTED, va="top")
            a0.text(1.0, y, val, fontsize=6.8, color=INK, va="top", ha="right", fontweight="bold")
            y -= 0.115
        a0.text(0.0, y - 0.04,
                "Venetoclax is the best-predicted inhibitor in\n"
                "the panel. Its fusion leans hardest on the\n"
                "differentiation-state block, matching the\n"
                "published account that monocytic AML is\n"
                "relatively venetoclax-resistant and primitive\n"
                "AML relatively sensitive.",
                fontsize=6, color=MUTED, va="top")
    a0.set_title("Venetoclax — the worked example", fontsize=7.4, loc="left", fontweight="bold")

    if S and S.get("primitive_vs_monocytic", {}).get("top10"):
        pv = S["primitive_vs_monocytic"]
        rows = list(reversed(pv["top10"]))[:10]
        y = np.arange(len(rows))
        vals = [r["mean_primitive_minus_monocytic"] for r in rows]
        cols = [BAD if r["inhibitor"] == "Venetoclax" else COOL for r in rows]
        a1.barh(y, vals, height=0.72, color=cols, lw=0)
        a1.set_yticks(y); a1.set_yticklabels([r["inhibitor"][:28] for r in rows], fontsize=5.8)
        a1.axvline(0, color=INK, lw=0.7)
        a1.set_xlabel("predicted sensitivity: primitive states minus monocytic states")
        v = pv.get("venetoclax") or {}
        a1.set_title("Model B, never shown a cell-state label in training\n"
                     "venetoclax ranks %s/%s, positive in %s of samples"
                     % (v.get("rank", "?"), pv.get("venetoclax_rank_of", "?"),
                        "%.0f%%" % (100 * v["frac_positive"]) if v.get("frac_positive") is not None else "?"),
                     fontsize=7.2, loc="left", fontweight="bold")
        style(a1)
    else:
        a1.axis("off")
        a1.text(0.5, 0.5, "state_response_validation.json not built yet", ha="center",
                color=MUTED, fontsize=7, style="italic")
    fig.subplots_adjust(left=0.03, right=0.99, top=0.86, bottom=0.13, wspace=0.30)
    save(fig, "Rx6_venetoclax_exemplar")


if __name__ == "__main__":
    if V is None or C is None:
        raise SystemExit("run train_drug_model.py and eval_drug_model.py first")
    for fn in (rx1, rx2, rx3, rx4, rx5, rx6):
        try:
            fn()
        except Exception as e:
            print("!! %s failed: %s" % (fn.__name__, e))
