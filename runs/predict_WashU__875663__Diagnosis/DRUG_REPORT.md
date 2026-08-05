# MOSAIC-Rx drug prioritisation

*Predicted **ex-vivo** sensitivity from the BeatAML2 functional screen. This is a prioritisation for trial matching or laboratory validation, **not** a treatment recommendation.*

- inhibitors modelled: **118**, reported 118, abstained 0
- genes shared with the model: 14237
- distance from the BeatAML training distribution: 53th percentile
- differentiation axis (primitive - mature): 1th percentile of BeatAML
- cell states scored: 20

## approved in AML  (6 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Midostaurin | 0.462 | 0.416 | 0.803 | 0.545 | — | 0.336 | 0 |
| 2 | Quizartinib (AC220) | 0.348 | 0.364 | 0.855 | 0.040 | — | 0.330 | 0 |
| 3 | Gilteritinib | 0.266 | 0.147 | 0.850 | 0.014 | — | 0.299 | 1 |
| 4 | Azacytidine | 0.233 | 0.329 | 0.646 | 0.022 | — | 0.107 | 2 |
| 5 | Cytarabine | 0.230 | 0.310 | 0.824 | 0.051 | — | 0.076 | 2 |
| 6 | Venetoclax | 0.193 | 0.002 | 0.966 | 0.004 | — | 0.138 | 2 |

## approved, other indication (off-label in AML)  (36 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Trametinib (GSK1120212) | 0.771 | 0.955 | 0.903 | 0.997 | — | 0.469 | 1 |
| 2 | Dasatinib | 0.767 | 0.970 | 0.931 | 0.994 | — | 0.373 | 1 |
| 3 | Selumetinib (AZD6244) | 0.763 | 0.946 | 0.886 | 0.992 | — | 0.469 | 1 |
| 4 | Tivozanib (AV-951) | 0.760 | 0.976 | 0.941 | 0.998 | — | 0.312 | 1 |
| 5 | Rapamycin | 0.745 | 0.959 | 0.922 | 0.996 | — | 0.302 | 2 |
| 6 | Idelalisib | 0.704 | 0.833 | 0.759 | 1.000 | — | 0.476 | 1 |
| 7 | Panobinostat | 0.681 | 0.827 | 0.860 | 0.982 | — | 0.366 | 1 |
| 8 | Lenvatinib | 0.668 | 0.824 | 0.861 | 0.970 | — | 0.302 | 1 |
| 9 | Axitinib (AG-013736) | 0.668 | 0.811 | 0.862 | 1.000 | — | 0.293 | 1 |
| 10 | Ponatinib (AP24534) | 0.666 | 0.795 | 0.952 | 0.944 | — | 0.347 | 1 |

## clinical-trial agent  (36 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Linifanib (ABT-869) | 0.700 | 0.918 | 0.903 | 0.995 | — | 0.434 | 1 |
| 2 | Foretinib (XL880) | 0.681 | 0.896 | 0.908 | 1.000 | — | 0.373 | 1 |
| 3 | CI-1040 (PD184352) | 0.659 | 0.852 | 0.819 | 0.968 | — | 0.469 | 0 |
| 4 | OTX-015 | 0.655 | 0.874 | 0.779 | 1.000 | — | 0.391 | 1 |
| 5 | GDC-0941 | 0.651 | 0.870 | 0.825 | 1.000 | — | 0.342 | 3 |
| 6 | Elesclomol | 0.637 | 0.806 | 0.800 | 1.000 | — | 0.458 | 1 |
| 7 | Motesanib (AMG-706) | 0.633 | 0.841 | 0.824 | 1.000 | — | 0.318 | 1 |
| 8 | INK-128 | 0.628 | 0.842 | 0.793 | 1.000 | — | 0.302 | 2 |
| 9 | Cediranib (AZD2171) | 0.610 | 0.825 | 0.748 | 0.998 | — | 0.293 | 1 |
| 10 | MK-2206 | 0.609 | 0.773 | 0.767 | 1.000 | — | 0.395 | 1 |

## research-only compound  (40 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | GW-2580 | 0.669 | 0.795 | 0.863 | 0.990 | — | 0.688 | 1 |
| 2 | PRT062607 | 0.636 | 0.743 | 0.804 | 1.000 | — | 0.682 | 0 |
| 3 | PP242 | 0.633 | 0.890 | 0.816 | 0.996 | — | 0.302 | 2 |
| 4 | JQ1 | 0.606 | 0.815 | 0.763 | 0.997 | — | 0.391 | 2 |
| 5 | Doramapimod (BIRB 796) | 0.601 | 0.817 | 0.841 | 0.994 | — | 0.515 | 1 |
| 6 | KI20227 | 0.581 | 0.821 | 0.808 | 1.000 | — | 0.412 | 3 |
| 7 | A-674563 | 0.568 | 0.657 | 0.777 | 0.942 | — | 0.623 | 2 |
| 8 | JAK Inhibitor I | 0.553 | 0.699 | 0.700 | 0.982 | — | 0.450 | 0 |
| 9 | TG100-115 | 0.535 | 0.691 | 0.689 | 0.990 | — | 0.380 | 1 |
| 10 | JNJ-28312141 | 0.508 | 0.665 | 0.842 | 0.921 | — | 0.498 | 1 |

## Combination hypotheses (complementary cell-state coverage)

*BeatAML2 measured single agents only; no synergy is claimed.*

| pair | pathways | coverage A | coverage B | union | gain |
|---|---|---|---|---|---|
| Azacytidine + CYT387 | epigenetic / JAK_STAT | 0.49 | 0.51 | 1.00 | +0.49 |
| CYT387 + Crizotinib (PF-2341066) | JAK_STAT / RTK | 0.51 | 0.55 | 1.00 | +0.45 |
| Crizotinib (PF-2341066) + Midostaurin | RTK / FLT3 | 0.55 | 0.66 | 1.00 | +0.34 |
| Azacytidine + Midostaurin | epigenetic / FLT3 | 0.49 | 0.66 | 0.94 | +0.28 |
