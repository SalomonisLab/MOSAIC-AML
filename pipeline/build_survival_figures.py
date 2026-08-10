#!/usr/bin/env python3
"""Journal-spec panels for the survival layer, generated from the validation JSONs.

  Sv1  what it adds over age + ELN 2017 -- the only comparison that matters
  Sv2  Kaplan-Meier by predicted risk tertile on the sealed hold-out, plus calibration
  Sv3  the honest limit: group timing is accurate, individual timing is not

  python build_survival_figures.py -> deliverables/figures/Sv{1,2,3}_*.pdf/.png
"""
import os, json
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
D = os.path.join(ROOT, "deliverables")
FIGD = os.path.join(D, "figures"); os.makedirs(FIGD, exist_ok=True)
MM = 1 / 25.4
INK, MUTED, GRID = "#1a1a19", "#6b6b63", "#d9d9d4"
ACC, GOOD, BAD, COOL = "#8a6a18", "#1a6e1a", "#a03828", "#9ec9e2"
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
                     "font.size": 7, "axes.labelsize": 7, "xtick.labelsize": 6.5,
                     "ytick.labelsize": 6.5, "legend.fontsize": 6.5, "axes.linewidth": 0.6,
                     "pdf.fonttype": 42, "savefig.dpi": 600})

C = json.load(open(os.path.join(D, "survival_model_card.json")))
T = json.load(open(os.path.join(D, "survival_time_validation.json")))


def style(ax, xg=False, yg=True):
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


LABEL = {"age": "age", "eln": "ELN 2017", "age_eln": "age + ELN", "clin": "clinical block",
         "rna": "RNA", "state": "cell state", "mut": "mutations",
         "molecular": "molecular only", "full": "full (molecular + clinical)"}


def sv1():
    order = ["eln", "state", "mut", "age", "rna", "molecular", "age_eln", "clin", "full"]
    fig, (a0, a1) = plt.subplots(1, 2, figsize=(180 * MM, 78 * MM),
                                 gridspec_kw={"width_ratios": [1, 1.05]})
    cv = [C["oof"][k]["c_index"] for k in order]
    ho = [C["holdout"][k]["c_index"] for k in order]
    y = np.arange(len(order))
    cols = [ACC if k == "full" else (GOOD if k == "age_eln" else COOL) for k in order]
    a0.barh(y - 0.19, [v - 0.5 for v in cv], 0.36, left=0.5, color=cols, lw=0, label="cross-validated")
    a0.barh(y + 0.19, [v - 0.5 for v in ho], 0.36, left=0.5, color=cols, lw=0, alpha=0.5,
            label="sealed hold-out")
    for i, (c, h) in enumerate(zip(cv, ho)):
        a0.text(c + 0.004, i - 0.19, "%.3f" % c, va="center", fontsize=5.6, color=INK)
        a0.text(h + 0.004, i + 0.19, "%.3f" % h, va="center", fontsize=5.6, color=MUTED)
    a0.axvline(C["oof"]["age_eln"]["c_index"], color=GOOD, lw=1.0, ls="--")
    a0.set_yticks(y); a0.set_yticklabels([LABEL[k] for k in order], fontsize=6.4)
    a0.set_xlim(0.5, 0.80); a0.set_xlabel("Harrell C-index  (0.5 = no discrimination)")
    a0.set_title("who dies first\n%d patients, %d deaths; dashed line = the age + ELN baseline"
                 % (C["cohort"]["final_patients"], C["cohort"]["events"]),
                 fontsize=7.2, loc="left", fontweight="bold")
    a0.legend(frameon=False, loc="lower right", fontsize=6)
    style(a0, xg=True, yg=False)

    inc = C["incremental"]
    keys = ["state", "mut", "rna", "molecular", "clin", "full"]
    y2 = np.arange(len(keys))
    d = [inc[k]["delta_c"] for k in keys]
    lo = [inc[k]["ci95"][0] for k in keys]; hi = [inc[k]["ci95"][1] for k in keys]
    sig = [inc[k]["ci95"][0] > 0 for k in keys]
    a1.errorbar(d, y2, xerr=[np.array(d) - np.array(lo), np.array(hi) - np.array(d)],
                fmt="o", ms=5, lw=1.2, capsize=2.5,
                color=INK, ecolor=MUTED, mfc="none", zorder=3)
    for i, (dd, s) in enumerate(zip(d, sig)):
        a1.plot([dd], [i], "o", ms=6, color=ACC if s else "#c9c6b8", zorder=4)
        a1.text(hi[i] + 0.006, i, "%+.3f%s" % (dd, "  p=%.3f" % inc[keys[i]]["p_gt0"] if s else ""),
                va="center", fontsize=5.8, color=ACC if s else MUTED)
    a1.axvline(0, color=INK, lw=0.8)
    a1.set_yticks(y2); a1.set_yticklabels([LABEL[k] for k in keys], fontsize=6.4)
    a1.set_xlabel("change in C-index versus age + ELN 2017  (95% bootstrap CI over patients)")
    a1.set_title("what the molecular data actually adds\nonly the full fusion clears zero",
                 fontsize=7.2, loc="left", fontweight="bold")
    style(a1, xg=True, yg=False)
    fig.subplots_adjust(left=0.115, right=0.94, top=0.85, bottom=0.135, wspace=0.42)
    save(fig, "Sv1_survival_discrimination")


def sv2():
    fig, (a0, a1) = plt.subplots(1, 2, figsize=(180 * MM, 76 * MM))
    kmd = C["risk_groups"]["km"]
    for name, col in (("low", GOOD), ("intermediate", ACC), ("high", BAD)):
        k = kmd.get(name)
        if not k:
            continue
        ts = [0.0] + k["times"]; ss = [1.0] + k["survival"]
        a0.step(ts, ss, where="post", color=col, lw=1.6,
                label="%s risk  (n=%d, %d deaths)" % (name, k["n"], k["events"]))
    lr = C["risk_groups"]["logrank_low_vs_high"]
    a0.set_xlabel("years from diagnosis"); a0.set_ylabel("survival")
    a0.set_ylim(0, 1.02); a0.set_xlim(0, 5)
    a0.set_title("sealed hold-out, split by predicted risk tertile\nlow vs high log-rank p = %.1e"
                 % lr["p"], fontsize=7.2, loc="left", fontweight="bold")
    a0.legend(frameon=False, loc="upper right", fontsize=6)
    style(a0)

    for j, (h, col) in enumerate((("1y", GOOD), ("2y", ACC), ("5y", BAD))):
        c = T["calibration"].get(h)
        if not c:
            continue
        p = [b["predicted"] for b in c["bins"]]; o = [b["observed_km"] for b in c["bins"]]
        a1.plot(p, o, "o-", color=col, lw=1.3, ms=4.5,
                label="%s  (mean gap %.3f)" % (h, c["mean_abs_calibration_gap"]))
    a1.plot([0, 1], [0, 1], color=GRID, lw=0.9, ls="--")
    a1.set_xlabel("predicted P(alive)"); a1.set_ylabel("observed (Kaplan-Meier)")
    a1.set_xlim(0, 1); a1.set_ylim(0, 1)
    a1.set_title("does a stated probability mean what it says?\nquartiles of predicted survival, hold-out",
                 fontsize=7.2, loc="left", fontweight="bold")
    a1.legend(frameon=False, loc="upper left", fontsize=6)
    style(a1)
    fig.subplots_adjust(left=0.075, right=0.985, top=0.85, bottom=0.14, wspace=0.26)
    save(fig, "Sv2_risk_groups_and_calibration")


def sv3():
    fig, (a0, a1) = plt.subplots(1, 2, figsize=(180 * MM, 74 * MM),
                                 gridspec_kw={"width_ratios": [1, 1]})
    tert = T["median_survival_by_risk_tertile"]
    names = [k for k in ("low", "intermediate", "high") if k in tert]
    x = np.arange(len(names)); w = 0.36
    pred = [tert[k]["predicted_median_years"] or np.nan for k in names]
    obs = [tert[k]["observed_km_median_years"] or np.nan for k in names]
    a0.bar(x - w / 2, pred, w, color=ACC, label="predicted median")
    a0.bar(x + w / 2, obs, w, color=COOL, label="observed (Kaplan-Meier)")
    for i, k in enumerate(names):
        v = tert[k]
        if v["observed_km_median_years"] is None:
            a0.text(i + w / 2, 0.15, "not reached\n(>half alive)", ha="center", fontsize=5.6, color=MUTED)
        else:
            a0.text(i, max(pred[i], obs[i]) + 0.18, "error %.2f y" % v["error_years"],
                    ha="center", fontsize=6, color=GOOD, fontweight="bold")
    a0.set_xticks(x); a0.set_xticklabels(names); a0.set_ylabel("median survival (years)")
    a0.set_ylim(0, max(v for v in pred if v == v) * 1.22)
    a0.set_title("group timing is accurate\nmedian survival per predicted risk tertile",
                 fontsize=7.2, loc="left", fontweight="bold")
    a0.legend(frameon=False, loc="upper right", fontsize=6)
    style(a0)

    sp = T["individual_spread_within_risk_group"]
    names2 = [k for k in ("low", "intermediate", "high") if k in sp]
    y = np.arange(len(names2))
    for i, k in enumerate(names2):
        v = sp[k]
        a1.plot(v["survival_years_p10_p90"], [i, i], lw=6, color=COOL, solid_capstyle="butt")
        a1.plot([v["median"]], [i], "|", ms=16, mew=2.2, color=BAD)
        a1.text(v["survival_years_p10_p90"][1] + 0.06, i,
                "median %.1f y  ·  n=%d" % (v["median"], v["n_deaths"]),
                va="center", fontsize=6, color=MUTED)
    a1.set_yticks(y); a1.set_yticklabels(names2)
    a1.set_xlabel("actual years survived (deaths in the hold-out)")
    pp = T.get("per_patient_point_estimate") or {}
    a1.set_title("individual timing is not\n10th–90th percentile spread inside one risk band; "
                 "per-patient MAE %.2f y" % pp.get("mae_years", float("nan")),
                 fontsize=7.2, loc="left", fontweight="bold")
    a1.set_xlim(0, max(v["survival_years_p10_p90"][1] for v in sp.values()) + 1.1)
    style(a1, xg=True, yg=False)
    fig.text(0.005, 0.008, "This is why the report shows a survival curve and a risk group rather than "
             "a single predicted number of months for one person.", fontsize=5.8, color=MUTED)
    fig.subplots_adjust(left=0.085, right=0.985, top=0.85, bottom=0.175, wspace=0.30)
    save(fig, "Sv3_group_versus_individual_timing")


if __name__ == "__main__":
    for fn in (sv1, sv2, sv3):
        try:
            fn()
        except Exception as ex:
            print("!! %s failed: %s" % (fn.__name__, ex))
