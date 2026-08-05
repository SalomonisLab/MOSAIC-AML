# MOSAIC-Rx drug prioritisation

*Predicted **ex-vivo** sensitivity from the BeatAML2 functional screen. This is a prioritisation for trial matching or laboratory validation, **not** a treatment recommendation.*

- inhibitors modelled: **118**, reported 118, abstained 0
- genes shared with the model: 14237
- distance from the BeatAML training distribution: 37th percentile
- differentiation axis (primitive - mature): 24th percentile of BeatAML
- cell states scored: 41

## approved in AML  (6 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Gilteritinib | 0.857 | 0.917 | 0.850 | 0.950 | 1.000 | 0.590 | 0 |
| 2 | Quizartinib (AC220) | 0.837 | 0.877 | 0.855 | 0.996 | 1.000 | 0.517 | 0 |
| 3 | Midostaurin | 0.836 | 0.892 | 0.803 | 0.999 | 1.000 | 0.501 | 0 |
| 4 | Cytarabine | 0.804 | 0.869 | 0.824 | 0.913 | 1.000 | 0.585 | 1 |
| 5 | Venetoclax | 0.448 | 0.026 | 0.966 | 0.433 | 1.000 | 0.333 | 1 |
| 6 | Azacytidine | 0.390 | 0.416 | 0.646 | 0.072 | 0.329 | 0.629 | 1 |

## approved, other indication (off-label in AML)  (36 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Ponatinib (AP24534) | 0.847 | 0.987 | 0.952 | 0.997 | 1.000 | 0.530 | 1 |
| 2 | Rapamycin | 0.839 | 0.936 | 0.922 | 0.970 | 0.961 | 0.730 | 0 |
| 3 | Sorafenib | 0.821 | 0.954 | 0.919 | 0.995 | 1.000 | 0.481 | 1 |
| 4 | Ibrutinib (PCI-32765) | 0.819 | 0.919 | 0.825 | 1.000 | 1.000 | 0.614 | 0 |
| 5 | Cabozantinib | 0.818 | 0.938 | 0.901 | 0.960 | 1.000 | 0.542 | 0 |
| 6 | Bosutinib (SKI-606) | 0.812 | 0.937 | 0.862 | 0.996 | 1.000 | 0.510 | 1 |
| 7 | Sunitinib | 0.810 | 0.931 | 0.891 | 0.992 | 1.000 | 0.496 | 1 |
| 8 | Idelalisib | 0.808 | 0.836 | 0.759 | 0.997 | 1.000 | 0.805 | 0 |
| 9 | Pazopanib (GW786034) | 0.801 | 0.943 | 0.885 | 0.999 | 1.000 | 0.404 | 1 |
| 10 | Dasatinib | 0.797 | 0.976 | 0.931 | 0.828 | 1.000 | 0.440 | 1 |

## clinical-trial agent  (36 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Entospletinib (GS-9973) | 0.828 | 0.965 | 0.872 | 1.000 | 1.000 | 0.790 | 0 |
| 2 | LY-333531 | 0.808 | 0.925 | 0.871 | 0.993 | 1.000 | 0.771 | 0 |
| 3 | Linifanib (ABT-869) | 0.777 | 0.946 | 0.903 | 0.986 | 1.000 | 0.508 | 1 |
| 4 | Saracatinib (AZD0530) | 0.774 | 0.940 | 0.872 | 0.998 | 1.000 | 0.517 | 1 |
| 5 | KW-2449 | 0.774 | 0.951 | 0.913 | 0.959 | 1.000 | 0.507 | 1 |
| 6 | GDC-0941 | 0.773 | 0.888 | 0.825 | 0.978 | 1.000 | 0.691 | 1 |
| 7 | 17-AAG (Tanespimycin) | 0.770 | 0.935 | 0.894 | 0.958 | 1.000 | 0.530 | 1 |
| 8 | Crenolanib | 0.764 | 0.923 | 0.879 | 0.993 | 1.000 | 0.496 | 1 |
| 9 | CI-1040 (PD184352) | 0.757 | 0.888 | 0.819 | 0.992 | 1.000 | 0.587 | 0 |
| 10 | RAF265 (CHIR-265) | 0.754 | 0.902 | 0.831 | 1.000 | 1.000 | 0.528 | 0 |

## research-only compound  (40 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | PRT062607 | 0.770 | 0.885 | 0.804 | 0.998 | 1.000 | 0.790 | 0 |
| 2 | NVP-TAE684 | 0.740 | 0.958 | 0.907 | 0.970 | 1.000 | 0.365 | 1 |
| 3 | GW-2580 | 0.738 | 0.945 | 0.863 | 0.831 | 1.000 | 0.553 | 0 |
| 4 | PI-103 | 0.734 | 0.827 | 0.757 | 0.998 | 1.000 | 0.726 | 0 |
| 5 | A-674563 | 0.725 | 0.817 | 0.777 | 0.959 | 1.000 | 0.724 | 0 |
| 6 | JNJ-28312141 | 0.718 | 0.910 | 0.842 | 1.000 | 1.000 | 0.571 | 1 |
| 7 | Doramapimod (BIRB 796) | 0.711 | 0.919 | 0.841 | 0.997 | 1.000 | 0.528 | 1 |
| 8 | SU11274 | 0.699 | 0.870 | 0.809 | 0.960 | 1.000 | 0.432 | 1 |
| 9 | TG100-115 | 0.691 | 0.716 | 0.689 | 0.996 | 1.000 | 0.806 | 0 |
| 10 | PD173955 | 0.689 | 0.864 | 0.794 | 0.826 | 1.000 | 0.510 | 0 |

## Combination hypotheses (complementary cell-state coverage)

*BeatAML2 measured single agents only; no synergy is claimed.*

| pair | pathways | coverage A | coverage B | union | gain |
|---|---|---|---|---|---|
| Bortezomib (Velcade) + Lenalidomide | proteostasis / epigenetic | 0.56 | 0.53 | 1.00 | +0.44 |
| Bortezomib (Velcade) + Palbociclib | proteostasis / cell_cycle | 0.56 | 0.53 | 1.00 | +0.44 |
| Entrectinib + Lenalidomide | RTK / epigenetic | 0.67 | 0.53 | 1.00 | +0.33 |
| Entrectinib + Palbociclib | RTK / cell_cycle | 0.67 | 0.53 | 1.00 | +0.33 |
| Palbociclib + Panobinostat | cell_cycle / epigenetic | 0.53 | 0.69 | 1.00 | +0.31 |
| Bortezomib (Velcade) + Neratinib (HKI-272) | proteostasis / RTK | 0.56 | 0.74 | 1.00 | +0.26 |
| Neratinib (HKI-272) + Panobinostat | RTK / epigenetic | 0.74 | 0.69 | 1.00 | +0.26 |
| Entrectinib + Panobinostat | RTK / epigenetic | 0.67 | 0.69 | 0.89 | +0.20 |
