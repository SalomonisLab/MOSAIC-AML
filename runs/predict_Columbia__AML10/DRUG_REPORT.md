# MOSAIC-Rx drug prioritisation

*Predicted **ex-vivo** sensitivity from the BeatAML2 functional screen. This is a prioritisation for trial matching or laboratory validation, **not** a treatment recommendation.*

- inhibitors modelled: **118**, reported 118, abstained 0
- genes shared with the model: 14237
- distance from the BeatAML training distribution: 47th percentile
- differentiation axis (primitive - mature): 12th percentile of BeatAML
- cell states scored: 44

## approved in AML  (6 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Cytarabine | 0.772 | 0.852 | 0.824 | 0.863 | 1.000 | 0.487 | 1 |
| 2 | Midostaurin | 0.718 | 0.806 | 0.803 | 0.764 | 0.842 | 0.498 | 0 |
| 3 | Quizartinib (AC220) | 0.596 | 0.565 | 0.855 | 0.577 | 0.845 | 0.498 | 0 |
| 4 | Gilteritinib | 0.566 | 0.458 | 0.850 | 0.596 | 0.842 | 0.459 | 0 |
| 5 | Azacytidine | 0.447 | 0.428 | 0.646 | 0.248 | 0.842 | 0.449 | 1 |
| 6 | Venetoclax | 0.246 | 0.006 | 0.966 | 0.054 | 0.362 | 0.178 | 1 |

## approved, other indication (off-label in AML)  (36 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Ponatinib (AP24534) | 0.794 | 0.978 | 0.952 | 0.969 | 0.878 | 0.462 | 1 |
| 2 | Ibrutinib (PCI-32765) | 0.721 | 0.837 | 0.825 | 0.942 | 0.845 | 0.537 | 0 |
| 3 | Rapamycin | 0.701 | 0.924 | 0.922 | 0.954 | 0.752 | 0.280 | 1 |
| 4 | Crizotinib (PF-2341066) | 0.674 | 0.758 | 0.802 | 0.959 | 0.845 | 0.434 | 0 |
| 5 | Idelalisib | 0.672 | 0.778 | 0.759 | 0.955 | 0.845 | 0.424 | 0 |
| 6 | CYT387 | 0.622 | 0.689 | 0.703 | 0.939 | 0.845 | 0.390 | 0 |
| 7 | Trametinib (GSK1120212) | 0.610 | 0.867 | 0.903 | 0.746 | 0.510 | 0.477 | 1 |
| 8 | Entrectinib | 0.601 | 0.577 | 0.771 | 0.735 | 1.000 | 0.490 | 1 |
| 9 | Afatinib (BIBW-2992) | 0.589 | 0.690 | 0.718 | 0.758 | 0.845 | 0.354 | 0 |
| 10 | Gefitinib | 0.549 | 0.605 | 0.794 | 0.529 | 0.845 | 0.506 | 0 |

## clinical-trial agent  (36 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Elesclomol | 0.646 | 0.792 | 0.800 | 0.946 | 0.638 | 0.790 | 0 |
| 2 | GDC-0941 | 0.629 | 0.806 | 0.825 | 0.940 | 0.811 | 0.365 | 2 |
| 3 | Alisertib (MLN8237) | 0.607 | 0.682 | 0.728 | 0.803 | 0.845 | 0.680 | 0 |
| 4 | Motesanib (AMG-706) | 0.607 | 0.788 | 0.824 | 0.945 | 0.687 | 0.474 | 0 |
| 5 | AZD1480 | 0.601 | 0.656 | 0.680 | 0.963 | 0.923 | 0.470 | 0 |
| 6 | INK-128 | 0.600 | 0.792 | 0.793 | 0.958 | 0.779 | 0.280 | 2 |
| 7 | Flavopiridol | 0.592 | 0.729 | 0.695 | 0.928 | 0.665 | 0.658 | 0 |
| 8 | SNS-032 (BMS-387032) | 0.591 | 0.745 | 0.801 | 0.952 | 0.687 | 0.494 | 0 |
| 9 | BEZ235 | 0.568 | 0.777 | 0.765 | 0.945 | 0.687 | 0.298 | 1 |
| 10 | Volasertib (BI-6727) | 0.554 | 0.620 | 0.668 | 0.675 | 0.845 | 0.706 | 0 |

## research-only compound  (40 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | A-674563 | 0.679 | 0.814 | 0.777 | 0.986 | 1.000 | 0.445 | 0 |
| 2 | JAK Inhibitor I | 0.638 | 0.738 | 0.700 | 0.986 | 1.000 | 0.421 | 0 |
| 3 | KI20227 | 0.628 | 0.845 | 0.808 | 0.966 | 0.842 | 0.498 | 2 |
| 4 | PI-103 | 0.614 | 0.798 | 0.757 | 0.971 | 0.878 | 0.298 | 1 |
| 5 | AGI-6780 | 0.612 | 0.554 | 0.773 | 0.954 | 0.923 | 0.845 | 0 |
| 6 | PD173955 | 0.609 | 0.766 | 0.794 | 0.939 | 0.845 | 0.403 | 0 |
| 7 | JQ1 | 0.609 | 0.787 | 0.763 | 0.959 | 0.845 | 0.365 | 1 |
| 8 | PP242 | 0.549 | 0.837 | 0.816 | 0.922 | 0.594 | 0.280 | 1 |
| 9 | Doramapimod (BIRB 796) | 0.548 | 0.864 | 0.841 | 0.892 | 0.446 | 0.652 | 1 |
| 10 | JNJ-28312141 | 0.532 | 0.770 | 0.842 | 0.719 | 0.746 | 0.488 | 2 |

## Combination hypotheses (complementary cell-state coverage)

*BeatAML2 measured single agents only; no synergy is claimed.*

| pair | pathways | coverage A | coverage B | union | gain |
|---|---|---|---|---|---|
| Palbociclib + Panobinostat | cell_cycle / epigenetic | 0.45 | 0.55 | 1.00 | +0.45 |
| Bosutinib (SKI-606) + Palbociclib | RTK / cell_cycle | 0.42 | 0.45 | 0.87 | +0.42 |
| Imatinib + Palbociclib | RTK / cell_cycle | 0.63 | 0.45 | 1.00 | +0.37 |
| Axitinib (AG-013736) + Palbociclib | RTK / cell_cycle | 0.66 | 0.45 | 1.00 | +0.34 |
| Axitinib (AG-013736) + Lenalidomide | RTK / epigenetic | 0.66 | 0.69 | 1.00 | +0.31 |
| Bosutinib (SKI-606) + Lenalidomide | RTK / epigenetic | 0.42 | 0.69 | 1.00 | +0.31 |
| Cabozantinib + Lenalidomide | RTK / epigenetic | 0.63 | 0.69 | 1.00 | +0.31 |
| Imatinib + Lenalidomide | RTK / epigenetic | 0.63 | 0.69 | 1.00 | +0.31 |
