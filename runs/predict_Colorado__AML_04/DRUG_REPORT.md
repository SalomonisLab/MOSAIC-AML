# MOSAIC-Rx drug prioritisation

*Predicted **ex-vivo** sensitivity from the BeatAML2 functional screen. This is a prioritisation for trial matching or laboratory validation, **not** a treatment recommendation.*

- inhibitors modelled: **118**, reported 118, abstained 0
- genes shared with the model: 14237
- distance from the BeatAML training distribution: 72th percentile
- differentiation axis (primitive - mature): 92th percentile of BeatAML
- cell states scored: 29

## approved in AML  (6 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Venetoclax | 0.801 | 0.989 | 0.966 | 0.956 | 1.000 | 0.715 | 2 |
| 2 | Quizartinib (AC220) | 0.794 | 0.900 | 0.855 | 0.992 | 1.000 | 0.733 | 0 |
| 3 | Cytarabine | 0.754 | 0.898 | 0.824 | 0.976 | 1.000 | 0.560 | 1 |
| 4 | Gilteritinib | 0.618 | 0.894 | 0.850 | 0.632 | 0.598 | 0.695 | 0 |
| 5 | Midostaurin | 0.359 | 0.599 | 0.803 | 0.389 | 0.302 | 0.604 | 1 |
| 6 | Azacytidine | 0.317 | 0.453 | 0.646 | 0.251 | 0.306 | 0.535 | 2 |

## approved, other indication (off-label in AML)  (36 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Ponatinib (AP24534) | 0.770 | 0.968 | 0.952 | 0.968 | 1.000 | 0.692 | 0 |
| 2 | Sorafenib | 0.751 | 0.961 | 0.919 | 0.966 | 1.000 | 0.622 | 0 |
| 3 | Cabozantinib | 0.747 | 0.948 | 0.901 | 0.975 | 1.000 | 0.621 | 0 |
| 4 | Sunitinib | 0.740 | 0.936 | 0.891 | 0.968 | 1.000 | 0.620 | 0 |
| 5 | Lenvatinib | 0.733 | 0.933 | 0.861 | 0.995 | 1.000 | 0.595 | 0 |
| 6 | Erlotinib | 0.725 | 0.917 | 0.834 | 0.992 | 1.000 | 0.595 | 0 |
| 7 | Pazopanib (GW786034) | 0.720 | 0.896 | 0.885 | 0.983 | 1.000 | 0.577 | 0 |
| 8 | Regorafenib (BAY 73-4506) | 0.716 | 0.908 | 0.834 | 0.980 | 1.000 | 0.565 | 0 |
| 9 | Ibrutinib (PCI-32765) | 0.703 | 0.834 | 0.825 | 0.987 | 1.000 | 0.670 | 0 |
| 10 | Vargetef | 0.701 | 0.877 | 0.827 | 0.973 | 1.000 | 0.571 | 0 |

## clinical-trial agent  (36 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | KW-2449 | 0.726 | 0.945 | 0.913 | 0.966 | 1.000 | 0.755 | 0 |
| 2 | Entospletinib (GS-9973) | 0.705 | 0.932 | 0.872 | 0.983 | 1.000 | 0.676 | 1 |
| 3 | Pelitinib (EKB-569) | 0.702 | 0.950 | 0.893 | 0.988 | 1.000 | 0.595 | 0 |
| 4 | LY-333531 | 0.698 | 0.918 | 0.871 | 0.968 | 1.000 | 0.700 | 1 |
| 5 | Dovitinib (CHIR-258) | 0.697 | 0.946 | 0.917 | 0.969 | 1.000 | 0.584 | 0 |
| 6 | Foretinib (XL880) | 0.689 | 0.921 | 0.908 | 0.981 | 1.000 | 0.575 | 0 |
| 7 | RAF265 (CHIR-265) | 0.645 | 0.879 | 0.831 | 0.966 | 1.000 | 0.494 | 1 |
| 8 | Saracatinib (AZD0530) | 0.634 | 0.760 | 0.872 | 0.991 | 1.000 | 0.660 | 0 |
| 9 | GDC-0941 | 0.630 | 0.797 | 0.825 | 0.995 | 1.000 | 0.583 | 2 |
| 10 | Canertinib (CI-1033) | 0.629 | 0.789 | 0.773 | 0.956 | 1.000 | 0.672 | 1 |

## research-only compound  (40 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | NVP-TAE684 | 0.692 | 0.963 | 0.907 | 0.985 | 1.000 | 0.595 | 0 |
| 2 | PRT062607 | 0.656 | 0.866 | 0.804 | 0.995 | 1.000 | 0.676 | 0 |
| 3 | Nutlin 3a | 0.637 | 0.828 | 0.786 | 0.964 | 1.000 | 0.719 | 3 |
| 4 | SU11274 | 0.636 | 0.871 | 0.809 | 0.995 | 1.000 | 0.547 | 0 |
| 5 | PHA-665752 | 0.636 | 0.876 | 0.798 | 0.986 | 1.000 | 0.547 | 0 |
| 6 | JNJ-28312141 | 0.620 | 0.833 | 0.842 | 0.970 | 1.000 | 0.738 | 1 |
| 7 | GSK-1838705A | 0.579 | 0.785 | 0.744 | 0.956 | 1.000 | 0.713 | 1 |
| 8 | NF-kB Activation Inhibitor | 0.566 | 0.941 | 0.865 | 0.961 | 0.981 | 0.180 | 3 |
| 9 | PD173955 | 0.565 | 0.643 | 0.794 | 0.973 | 0.994 | 0.709 | 0 |
| 10 | AGI-6780 | 0.558 | 0.783 | 0.773 | 0.983 | 1.000 | 0.354 | 1 |

## Combination hypotheses (complementary cell-state coverage)

*BeatAML2 measured single agents only; no synergy is claimed.*

| pair | pathways | coverage A | coverage B | union | gain |
|---|---|---|---|---|---|
| Axitinib (AG-013736) + Midostaurin | RTK / FLT3 | 0.57 | 0.65 | 1.00 | +0.35 |
