# MOSAIC-Rx drug prioritisation

*Predicted **ex-vivo** sensitivity from the BeatAML2 functional screen. This is a prioritisation for trial matching or laboratory validation, **not** a treatment recommendation.*

- inhibitors modelled: **118**, reported 118, abstained 0
- genes shared with the model: 14237
- distance from the BeatAML training distribution: 90th percentile
- differentiation axis (primitive - mature): 0th percentile of BeatAML
- cell states scored: 23

## approved in AML  (6 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Gilteritinib | 0.221 | 0.506 | 0.850 | 0.667 | 0.000 | 0.595 | 1 |
| 2 | Quizartinib (AC220) | 0.153 | 0.325 | 0.855 | 0.071 | 0.000 | 0.505 | 2 |
| 3 | Midostaurin | 0.150 | 0.296 | 0.803 | 0.158 | 0.000 | 0.481 | 1 |
| 4 | Venetoclax | 0.135 | 0.002 | 0.966 | 0.003 | 1.000 | 0.093 | 3 |
| 5 | Azacytidine | 0.010 | 0.423 | 0.646 | 0.066 | 0.000 | 0.000 | 3 |
| 6 | Cytarabine | -0.058 | 0.150 | 0.824 | 0.008 | 0.000 | 0.000 | 2 |

## approved, other indication (off-label in AML)  (36 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Selumetinib (AZD6244) | 0.682 | 0.948 | 0.886 | 0.979 | 1.000 | 0.705 | 1 |
| 2 | Gefitinib | 0.629 | 0.752 | 0.794 | 0.991 | 1.000 | 0.884 | 1 |
| 3 | Ibrutinib (PCI-32765) | 0.623 | 0.916 | 0.825 | 0.982 | 1.000 | 0.456 | 2 |
| 4 | Panobinostat | 0.552 | 0.857 | 0.860 | 0.991 | 1.000 | 0.140 | 2 |
| 5 | Nilotinib | 0.540 | 0.675 | 0.690 | 0.963 | 1.000 | 0.647 | 2 |
| 6 | Lapatinib | 0.505 | 0.697 | 0.734 | 0.871 | 1.000 | 0.451 | 2 |
| 7 | Tofacitinib (CP-690550) | 0.487 | 0.682 | 0.660 | 0.964 | 1.000 | 0.339 | 1 |
| 8 | CYT387 | 0.481 | 0.630 | 0.703 | 0.819 | 1.000 | 0.540 | 1 |
| 9 | Ruxolitinib (INCB018424) | 0.455 | 0.563 | 0.713 | 0.819 | 1.000 | 0.550 | 1 |
| 10 | Trametinib (GSK1120212) | 0.423 | 0.941 | 0.903 | 0.970 | 0.000 | 0.705 | 2 |

## clinical-trial agent  (36 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | CI-1040 (PD184352) | 0.614 | 0.889 | 0.819 | 0.986 | 1.000 | 0.705 | 1 |
| 2 | Saracatinib (AZD0530) | 0.585 | 0.801 | 0.872 | 0.953 | 1.000 | 0.763 | 1 |
| 3 | GDC-0941 | 0.559 | 0.899 | 0.825 | 0.992 | 1.000 | 0.352 | 4 |
| 4 | BEZ235 | 0.524 | 0.829 | 0.765 | 0.993 | 1.000 | 0.339 | 3 |
| 5 | MK-2206 | 0.400 | 0.784 | 0.767 | 0.992 | 1.000 | 0.258 | 3 |
| 6 | Volasertib (BI-6727) | 0.385 | 0.649 | 0.668 | 0.926 | 1.000 | 0.092 | 3 |
| 7 | Linifanib (ABT-869) | 0.362 | 0.895 | 0.903 | 0.946 | 0.000 | 0.706 | 2 |
| 8 | Flavopiridol | 0.358 | 0.778 | 0.695 | 0.993 | 1.000 | 0.062 | 3 |
| 9 | Foretinib (XL880) | 0.331 | 0.948 | 0.908 | 0.997 | 0.000 | 0.681 | 2 |
| 10 | Motesanib (AMG-706) | 0.323 | 0.825 | 0.824 | 0.970 | 0.000 | 0.683 | 1 |

## research-only compound  (40 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | PP242 | 0.528 | 0.900 | 0.816 | 0.991 | 1.000 | 0.267 | 3 |
| 2 | Doramapimod (BIRB 796) | 0.519 | 0.846 | 0.841 | 0.963 | 1.000 | 0.574 | 2 |
| 3 | PI-103 | 0.510 | 0.839 | 0.757 | 0.999 | 1.000 | 0.339 | 3 |
| 4 | TG100-115 | 0.464 | 0.688 | 0.689 | 0.990 | 1.000 | 0.501 | 2 |
| 5 | H-89 | 0.454 | 0.600 | 0.667 | 0.958 | 1.000 | 0.891 | 2 |
| 6 | STO609 | 0.320 | 0.493 | 0.612 | 0.607 | 1.000 | 0.732 | 3 |
| 7 | JNJ-7706621 | 0.292 | 0.531 | 0.663 | 0.736 | 1.000 | 0.105 | 3 |
| 8 | PD173955 | 0.289 | 0.785 | 0.794 | 0.969 | 0.000 | 0.691 | 2 |
| 9 | KI20227 | 0.284 | 0.873 | 0.808 | 0.970 | 0.000 | 0.646 | 3 |
| 10 | A-674563 | 0.269 | 0.788 | 0.777 | 0.980 | 0.000 | 0.564 | 3 |
