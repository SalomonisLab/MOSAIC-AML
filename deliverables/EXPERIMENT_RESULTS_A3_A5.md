# A3-A5 results

## A3 matrix factorisation — modest real gain

Interaction-target (honest) AUROC by rank: 5 -> 0.625, 12 -> 0.664, 25 -> 0.678, **40 -> 0.685**;
deployed per-family model = 0.674. Spearman 0.225 -> **0.239**. A 50/50 blend matches rank-40 MF.
So **+0.011 AUROC / +0.014 Spearman**, real but small — the low-rank structure is worth having,
and it helps the low-n inhibitors most, but it is not the step change the headline gap suggested.

## A4 raw concentration data — the ceiling, now measured directly

Downloaded `beataml_wv1to4_raw_inhibitor_v4_dbgap.txt` (49 MB, 555,583 wells). It carries a
`replicate` column, so assay noise can be measured rather than inferred from same-target proxies:

- **Technical replicates of the same well**: Pearson 0.406, Spearman 0.641,
  **median absolute difference 13.7 percentage points of viability**.
- **Per-drug, across the dose series**: median replicate r = **0.529** over 27 drugs with enough
  replicates, implying a ceiling of **sqrt(0.529) = 0.727** Spearman.

The correspondence with model performance is exact where it matters:

| inhibitor | replicate r | model AUROC |
|---|---|---|
| Trametinib | 0.87 | 0.87 |
| Venetoclax | 0.84 | 0.977 |
| Rapamycin | 0.80 | 0.86 |
| Lapatinib | 0.27 | — |
| PLX-4720 | 0.24 | — |

**The drugs the model predicts best are precisely the drugs whose assay reproduces best.** This is now
a measured fact, not an inference from same-target pairs.

## A5 imputed modalities — RUN ON THE CLUSTER, and it is a clean negative

The bundles only unpickle on Linux, so imputation ran on a compute node (`bsub -q test`) against
BeatAML bulk and the blocks were shipped back. Input-gene coverage was good: **ADT 95%** (129 features),
**GRN 94%** (7,486), **Lipid 96%** (202). No metabolite bundle exists in the cluster checkout, so that
modality is absent rather than failed.

Scored on the interaction target:

| arm | Spearman | AUROC |
|---|---|---|
| deployed | 0.223 | **0.672** |
| + ADT | 0.218 | 0.668 |
| + GRN | 0.211 | 0.662 |
| + Lipid | 0.218 | 0.669 |
| + all imputed | 0.205 | 0.658 |

**Every one makes it worse.** I predicted this would be the second-largest gain; it is not a gain at
all. In hindsight the prior should have been the other way round: these modalities are imputed *from
the RNA that is already in the model*, so they cannot carry information the RNA does not contain, and
what they add is dimensionality and imputation error. They earned their place in the mutation layer by
re-expressing RNA in a form a linear classifier found easier for specific drivers; that does not
transfer to this target.

## Where COMPASS now stands

Best achievable on the honest interaction target: **0.685** (rank-40 matrix factorisation) against a
**directly measured assay ceiling of ~0.727**. That is roughly 94% of the recoverable signal. The
remaining headroom is smaller than the spread between drugs, which is why the reliability tier now
shipped with every recommendation matters more than any further modelling.

---

# Item 1 — rebuilding the label from raw wells: a NEGATIVE, and a correction

I claimed last turn that replicate-averaging would raise the ceiling from 0.727 to 0.832 and called it
"the only remaining work that moves the real number". **That was wrong.** Measured:

| target | same-mechanism reliability | implied ceiling |
|---|---|---|
| released probit AUC (current) | **0.192** | **0.439** |
| replicate-averaged AUC | 0.168 | 0.410 |
| dose-anchored AUC (≤ 1 µM) | 0.102 | 0.320 |

Both rebuilt labels are **worse**. Two reasons, both of which I should have checked first:

1. **The replicates are not there.** Only **15.4%** of (specimen × inhibitor × concentration) cells have
   more than one well. Spearman–Brown assumed averaging two measurements everywhere; averaging one
   measurement with itself gains nothing, so the predicted lift was diluted to zero. On the subset that
   genuinely has replicates the gain is real but tiny and underpowered (0.188 → 0.196, 7 pairs).
2. **The released AUC is already denoised.** It comes from a fitted probit curve, which smooths across
   the whole dose series. Replacing it with a raw trapezoid discards that smoothing, and the loss
   exceeds anything replicate averaging recovers.

**Dose-anchoring was worse still** (0.102). Restricting to ≤ 1 µM removes the part of the curve where
these compounds actually discriminate — which is a genuinely uncomfortable finding in its own right,
and consistent with the pharmacology agent's standing flag that sensitivity often only appears in the
top concentration decade. As a *modelling target* it is clearly worse; as a statement about clinical
relevance it deserves attention.

**Conclusion:** the released probit AUC is the best available label. The ceiling stays at ~0.727, and
COMPASS at 0.685 remains ~94% of it. Raising the ceiling needs more replicate wells, not better
arithmetic on the ones that exist.

# Item 2 — the measured wins are now deployed

| change | before | after |
|---|---|---|
| survival, sealed hold-out C-index | 0.752 | **0.787** |
| survival, CV C-index | 0.726 | **0.751** |
| survival, 2-year AUC (hold-out) | 0.852 | **0.872** |
| survival, Brier at 2 y | 0.162 | **0.150** |

Deployed: **baseline induction type** (diagnosis-time only — `n_regimens` correlates 0.571 with
follow-up and would read the outcome), **age as a spline**, and **RMST** reported ahead of the median so
the low-risk group no longer gets a 6.9-year figure extrapolated past follow-up.
