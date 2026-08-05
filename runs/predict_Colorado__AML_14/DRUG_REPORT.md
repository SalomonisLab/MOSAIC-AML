# MOSAIC-Rx drug prioritisation

*Predicted **ex-vivo** sensitivity from the BeatAML2 functional screen. This is a prioritisation for trial matching or laboratory validation, **not** a treatment recommendation.*

- inhibitors modelled: **118**, reported 118, abstained 0
- genes shared with the model: 14237
- distance from the BeatAML training distribution: 57th percentile
- differentiation axis (primitive - mature): 53th percentile of BeatAML
- cell states scored: 37

## approved in AML  (6 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Midostaurin | 0.810 | 0.859 | 0.803 | 0.966 | 1.000 | 0.639 | 0 |
| 2 | Quizartinib (AC220) | 0.775 | 0.753 | 0.855 | 0.980 | 0.980 | 0.656 | 0 |
| 3 | Gilteritinib | 0.730 | 0.710 | 0.850 | 0.761 | 0.980 | 0.738 | 1 |
| 4 | Cytarabine | 0.579 | 0.649 | 0.824 | 0.648 | 0.709 | 0.682 | 1 |
| 5 | Venetoclax | 0.480 | 0.034 | 0.966 | 0.587 | 1.000 | 0.726 | 2 |
| 6 | Azacytidine | 0.302 | 0.301 | 0.646 | 0.000 | 0.000 | 0.857 | 3 |

## approved, other indication (off-label in AML)  (36 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Tivozanib (AV-951) | 0.844 | 0.973 | 0.941 | 0.970 | 1.000 | 0.747 | 0 |
| 2 | Dasatinib | 0.805 | 0.939 | 0.931 | 0.972 | 1.000 | 0.610 | 0 |
| 3 | Axitinib (AG-013736) | 0.787 | 0.867 | 0.862 | 1.000 | 1.000 | 0.688 | 0 |
| 4 | Pazopanib (GW786034) | 0.776 | 0.857 | 0.885 | 0.984 | 1.000 | 0.646 | 0 |
| 5 | Sorafenib | 0.755 | 0.875 | 0.919 | 0.796 | 1.000 | 0.661 | 0 |
| 6 | Entrectinib | 0.743 | 0.825 | 0.771 | 1.000 | 1.000 | 0.619 | 1 |
| 7 | Lenvatinib | 0.737 | 0.874 | 0.861 | 0.903 | 0.892 | 0.650 | 0 |
| 8 | Cabozantinib | 0.728 | 0.791 | 0.901 | 0.900 | 0.980 | 0.634 | 1 |
| 9 | Ponatinib (AP24534) | 0.694 | 0.774 | 0.952 | 0.688 | 0.980 | 0.662 | 1 |
| 10 | Panobinostat | 0.690 | 0.842 | 0.860 | 0.872 | 0.790 | 0.664 | 0 |

## clinical-trial agent  (36 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | OTX-015 | 0.728 | 0.856 | 0.779 | 0.984 | 1.000 | 0.708 | 0 |
| 2 | Saracatinib (AZD0530) | 0.717 | 0.834 | 0.872 | 0.954 | 1.000 | 0.654 | 0 |
| 3 | Entospletinib (GS-9973) | 0.705 | 0.892 | 0.872 | 0.885 | 0.892 | 0.698 | 0 |
| 4 | LY-333531 | 0.688 | 0.748 | 0.871 | 0.862 | 1.000 | 0.808 | 1 |
| 5 | AT7519 | 0.682 | 0.759 | 0.787 | 0.984 | 1.000 | 0.657 | 0 |
| 6 | Linifanib (ABT-869) | 0.663 | 0.811 | 0.903 | 0.857 | 0.872 | 0.663 | 0 |
| 7 | Crenolanib | 0.629 | 0.786 | 0.879 | 0.630 | 1.000 | 0.573 | 0 |
| 8 | Flavopiridol | 0.589 | 0.673 | 0.695 | 0.910 | 0.872 | 0.637 | 0 |
| 9 | AZD1480 | 0.518 | 0.589 | 0.680 | 0.835 | 0.790 | 0.661 | 0 |
| 10 | YM-155 | 0.515 | 0.723 | 0.754 | 0.911 | 0.808 | 0.671 | 0 |

## research-only compound  (40 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | PRT062607 | 0.693 | 0.821 | 0.804 | 0.968 | 1.000 | 0.698 | 0 |
| 2 | JNJ-28312141 | 0.667 | 0.822 | 0.842 | 0.936 | 1.000 | 0.721 | 1 |
| 3 | KI20227 | 0.647 | 0.811 | 0.808 | 0.946 | 0.980 | 0.670 | 1 |
| 4 | A-674563 | 0.566 | 0.648 | 0.777 | 0.815 | 0.980 | 0.545 | 1 |
| 5 | Nutlin 3a | 0.550 | 0.636 | 0.786 | 0.748 | 0.872 | 0.737 | 2 |
| 6 | GSK-1838705A | 0.422 | 0.447 | 0.744 | 0.560 | 1.000 | 0.590 | 2 |
| 7 | JNJ-7706621 | 0.418 | 0.504 | 0.663 | 0.679 | 0.626 | 0.767 | 0 |
| 8 | GW-2580 | 0.331 | 0.460 | 0.863 | 0.453 | 0.102 | 0.734 | 1 |
| 9 | NF-kB Activation Inhibitor | 0.318 | 0.322 | 0.865 | 0.527 | 0.851 | 0.248 | 2 |
| 10 | JQ1 | 0.305 | 0.576 | 0.763 | 0.558 | 0.149 | 0.708 | 1 |

## Combination hypotheses (complementary cell-state coverage)

*BeatAML2 measured single agents only; no synergy is claimed.*

| pair | pathways | coverage A | coverage B | union | gain |
|---|---|---|---|---|---|
| Azacytidine + Bortezomib (Velcade) | epigenetic / proteostasis | 0.63 | 0.55 | 1.00 | +0.37 |
| Azacytidine + CYT387 | epigenetic / JAK_STAT | 0.63 | 0.60 | 1.00 | +0.37 |
| Azacytidine + Ibrutinib (PCI-32765) | epigenetic / immune_signalling | 0.63 | 0.45 | 0.98 | +0.35 |
| Azacytidine + Idelalisib | epigenetic / PI3K_AKT_mTOR | 0.63 | 0.42 | 0.98 | +0.35 |
| Bortezomib (Velcade) + Palbociclib | proteostasis / cell_cycle | 0.55 | 0.67 | 1.00 | +0.33 |
| CYT387 + Palbociclib | JAK_STAT / cell_cycle | 0.60 | 0.67 | 1.00 | +0.33 |
| Ibrutinib (PCI-32765) + Lapatinib | immune_signalling / RTK | 0.45 | 0.53 | 0.85 | +0.33 |
| Ibrutinib (PCI-32765) + Palbociclib | immune_signalling / cell_cycle | 0.45 | 0.67 | 1.00 | +0.33 |
