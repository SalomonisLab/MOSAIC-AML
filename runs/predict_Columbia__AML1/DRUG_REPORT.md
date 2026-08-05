# MOSAIC-Rx drug prioritisation

*Predicted **ex-vivo** sensitivity from the BeatAML2 functional screen. This is a prioritisation for trial matching or laboratory validation, **not** a treatment recommendation.*

- inhibitors modelled: **118**, reported 118, abstained 0
- genes shared with the model: 14237
- distance from the BeatAML training distribution: 12th percentile
- differentiation axis (primitive - mature): 33th percentile of BeatAML
- cell states scored: 45

## approved in AML  (6 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Gilteritinib | 0.845 | 0.876 | 0.850 | 0.992 | 1.000 | 0.572 | 0 |
| 2 | Quizartinib (AC220) | 0.820 | 0.823 | 0.855 | 0.979 | 1.000 | 0.564 | 0 |
| 3 | Midostaurin | 0.549 | 0.670 | 0.803 | 0.501 | 0.565 | 0.501 | 0 |
| 4 | Cytarabine | 0.521 | 0.652 | 0.824 | 0.545 | 0.785 | 0.118 | 2 |
| 5 | Venetoclax | 0.375 | 0.015 | 0.966 | 0.375 | 0.785 | 0.184 | 1 |
| 6 | Azacytidine | 0.192 | 0.411 | 0.646 | 0.001 | 0.000 | 0.028 | 2 |

## approved, other indication (off-label in AML)  (36 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Rapamycin | 0.815 | 0.938 | 0.922 | 0.989 | 1.000 | 0.518 | 1 |
| 2 | Ponatinib (AP24534) | 0.788 | 0.978 | 0.952 | 0.906 | 0.819 | 0.583 | 0 |
| 3 | Trametinib (GSK1120212) | 0.771 | 0.852 | 0.903 | 0.980 | 0.986 | 0.494 | 0 |
| 4 | Vargetef | 0.754 | 0.778 | 0.827 | 0.913 | 1.000 | 0.676 | 0 |
| 5 | Lenvatinib | 0.754 | 0.897 | 0.861 | 0.900 | 0.802 | 0.673 | 0 |
| 6 | Selumetinib (AZD6244) | 0.749 | 0.815 | 0.886 | 0.964 | 0.982 | 0.494 | 0 |
| 7 | Regorafenib (BAY 73-4506) | 0.735 | 0.875 | 0.834 | 0.895 | 0.802 | 0.650 | 0 |
| 8 | Dasatinib | 0.734 | 0.901 | 0.931 | 0.886 | 0.785 | 0.543 | 0 |
| 9 | Bosutinib (SKI-606) | 0.716 | 0.826 | 0.862 | 0.890 | 0.785 | 0.644 | 0 |
| 10 | Cabozantinib | 0.716 | 0.878 | 0.901 | 0.979 | 1.000 | 0.622 | 0 |

## clinical-trial agent  (36 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Foretinib (XL880) | 0.796 | 0.907 | 0.908 | 0.985 | 1.000 | 0.719 | 0 |
| 2 | Entospletinib (GS-9973) | 0.738 | 0.916 | 0.872 | 0.913 | 1.000 | 0.464 | 1 |
| 3 | INK-128 | 0.723 | 0.841 | 0.793 | 0.994 | 1.000 | 0.518 | 0 |
| 4 | GDC-0941 | 0.711 | 0.847 | 0.825 | 0.979 | 0.986 | 0.453 | 2 |
| 5 | Saracatinib (AZD0530) | 0.692 | 0.863 | 0.872 | 0.895 | 0.802 | 0.648 | 0 |
| 6 | 17-AAG (Tanespimycin) | 0.670 | 0.915 | 0.894 | 0.883 | 0.777 | 0.420 | 1 |
| 7 | Motesanib (AMG-706) | 0.654 | 0.776 | 0.824 | 0.895 | 0.802 | 0.659 | 0 |
| 8 | Linifanib (ABT-869) | 0.641 | 0.798 | 0.903 | 0.811 | 0.785 | 0.579 | 0 |
| 9 | Cediranib (AZD2171) | 0.623 | 0.744 | 0.748 | 0.888 | 0.785 | 0.638 | 0 |
| 10 | BEZ235 | 0.617 | 0.806 | 0.765 | 0.891 | 0.763 | 0.454 | 0 |

## research-only compound  (40 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | PI-103 | 0.681 | 0.827 | 0.757 | 0.985 | 1.000 | 0.454 | 0 |
| 2 | KI20227 | 0.663 | 0.798 | 0.808 | 0.985 | 1.000 | 0.572 | 1 |
| 3 | A-674563 | 0.624 | 0.690 | 0.777 | 0.988 | 1.000 | 0.402 | 0 |
| 4 | JNJ-28312141 | 0.586 | 0.824 | 0.842 | 0.879 | 0.789 | 0.450 | 1 |
| 5 | PD173955 | 0.581 | 0.720 | 0.794 | 0.746 | 0.785 | 0.644 | 0 |
| 6 | PP242 | 0.564 | 0.845 | 0.816 | 0.787 | 0.565 | 0.518 | 0 |
| 7 | Doramapimod (BIRB 796) | 0.559 | 0.719 | 0.841 | 0.876 | 0.761 | 0.609 | 1 |
| 8 | JNJ-38877605 | 0.506 | 0.674 | 0.709 | 0.904 | 0.802 | 0.713 | 0 |
| 9 | AGI-5198 | 0.492 | 0.536 | 0.649 | 0.884 | 0.781 | 0.544 | 1 |
| 10 | JAK Inhibitor I | 0.487 | 0.684 | 0.700 | 0.814 | 0.615 | 0.421 | 0 |

## Combination hypotheses (complementary cell-state coverage)

*BeatAML2 measured single agents only; no synergy is claimed.*

| pair | pathways | coverage A | coverage B | union | gain |
|---|---|---|---|---|---|
| Bortezomib (Velcade) + Palbociclib | proteostasis / cell_cycle | 0.53 | 0.51 | 0.90 | +0.37 |
| CYT387 + Palbociclib | JAK_STAT / cell_cycle | 0.53 | 0.51 | 0.90 | +0.37 |
| Palbociclib + Panobinostat | cell_cycle / epigenetic | 0.51 | 0.53 | 0.90 | +0.37 |
| Bortezomib (Velcade) + Erlotinib | proteostasis / RTK | 0.53 | 0.78 | 1.00 | +0.22 |
| CYT387 + Erlotinib | JAK_STAT / RTK | 0.53 | 0.78 | 1.00 | +0.22 |
| Erlotinib + Panobinostat | RTK / epigenetic | 0.78 | 0.53 | 1.00 | +0.22 |
