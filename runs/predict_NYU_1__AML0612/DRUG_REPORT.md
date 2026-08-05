# MOSAIC-Rx drug prioritisation

*Predicted **ex-vivo** sensitivity from the BeatAML2 functional screen. This is a prioritisation for trial matching or laboratory validation, **not** a treatment recommendation.*

- inhibitors modelled: **118**, reported 118, abstained 0
- genes shared with the model: 14237
- distance from the BeatAML training distribution: 38th percentile
- differentiation axis (primitive - mature): 96th percentile of BeatAML
- cell states scored: 24

## approved in AML  (6 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Venetoclax | 0.907 | 0.973 | 0.966 | 0.981 | 1.000 | 0.632 | 1 |
| 2 | Quizartinib (AC220) | 0.805 | 0.832 | 0.855 | 0.976 | 1.000 | 0.442 | 0 |
| 3 | Gilteritinib | 0.663 | 0.645 | 0.850 | 0.847 | 0.859 | 0.420 | 0 |
| 4 | Midostaurin | 0.389 | 0.337 | 0.803 | 0.302 | 0.297 | 0.437 | 0 |
| 5 | Cytarabine | 0.288 | 0.398 | 0.824 | 0.170 | 0.163 | 0.264 | 1 |
| 6 | Azacytidine | 0.152 | 0.368 | 0.646 | 0.001 | 0.001 | 0.028 | 1 |

## approved, other indication (off-label in AML)  (36 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Erlotinib | 0.773 | 0.895 | 0.834 | 1.000 | 1.000 | 0.378 | 1 |
| 2 | Regorafenib (BAY 73-4506) | 0.768 | 0.863 | 0.834 | 0.976 | 0.987 | 0.473 | 0 |
| 3 | Crizotinib (PF-2341066) | 0.767 | 0.882 | 0.802 | 0.989 | 1.000 | 0.423 | 1 |
| 4 | Vargetef | 0.745 | 0.800 | 0.827 | 0.986 | 1.000 | 0.462 | 0 |
| 5 | Gefitinib | 0.715 | 0.772 | 0.794 | 0.987 | 1.000 | 0.378 | 0 |
| 6 | Palbociclib | 0.714 | 0.781 | 0.782 | 0.990 | 1.000 | 0.374 | 2 |
| 7 | Lapatinib | 0.711 | 0.793 | 0.734 | 0.986 | 1.000 | 0.352 | 1 |
| 8 | Lenvatinib | 0.675 | 0.653 | 0.861 | 0.963 | 0.987 | 0.468 | 0 |
| 9 | Entrectinib | 0.667 | 0.705 | 0.771 | 0.975 | 1.000 | 0.345 | 1 |
| 10 | Ruxolitinib (INCB018424) | 0.667 | 0.634 | 0.713 | 0.987 | 1.000 | 0.493 | 0 |

## clinical-trial agent  (36 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | RAF265 (CHIR-265) | 0.748 | 0.848 | 0.831 | 1.000 | 1.000 | 0.598 | 0 |
| 2 | Pelitinib (EKB-569) | 0.739 | 0.906 | 0.893 | 0.983 | 1.000 | 0.378 | 1 |
| 3 | Saracatinib (AZD0530) | 0.723 | 0.853 | 0.872 | 0.977 | 1.000 | 0.441 | 1 |
| 4 | Canertinib (CI-1033) | 0.676 | 0.800 | 0.773 | 0.980 | 1.000 | 0.352 | 1 |
| 5 | Masitinib (AB-1010) | 0.655 | 0.720 | 0.684 | 0.989 | 1.000 | 0.466 | 1 |
| 6 | Entospletinib (GS-9973) | 0.643 | 0.688 | 0.872 | 0.962 | 0.987 | 0.409 | 0 |
| 7 | Tozasertib (VX-680) | 0.635 | 0.700 | 0.728 | 1.000 | 1.000 | 0.360 | 1 |
| 8 | SNS-032 (BMS-387032) | 0.616 | 0.681 | 0.801 | 0.976 | 0.987 | 0.302 | 2 |
| 9 | GDC-0941 | 0.612 | 0.609 | 0.825 | 0.974 | 0.987 | 0.447 | 3 |
| 10 | Volasertib (BI-6727) | 0.588 | 0.684 | 0.668 | 0.990 | 1.000 | 0.192 | 2 |

## research-only compound  (40 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | NVP-TAE684 | 0.740 | 0.937 | 0.907 | 1.000 | 1.000 | 0.378 | 1 |
| 2 | NF-kB Activation Inhibitor | 0.727 | 0.921 | 0.865 | 0.984 | 1.000 | 0.590 | 1 |
| 3 | GW-2580 | 0.703 | 0.925 | 0.863 | 0.989 | 1.000 | 0.234 | 1 |
| 4 | PHA-665752 | 0.697 | 0.854 | 0.798 | 0.983 | 1.000 | 0.423 | 1 |
| 5 | PD173955 | 0.690 | 0.810 | 0.794 | 0.989 | 1.000 | 0.477 | 0 |
| 6 | ABT-737 | 0.677 | 0.776 | 0.763 | 0.960 | 0.987 | 0.596 | 2 |
| 7 | Nutlin 3a | 0.665 | 0.727 | 0.786 | 0.968 | 0.987 | 0.593 | 2 |
| 8 | JNJ-28312141 | 0.660 | 0.830 | 0.842 | 0.988 | 1.000 | 0.429 | 2 |
| 9 | SU11274 | 0.635 | 0.703 | 0.809 | 0.984 | 1.000 | 0.423 | 0 |
| 10 | GSK-1838705A | 0.627 | 0.801 | 0.744 | 0.999 | 1.000 | 0.369 | 2 |

## Combination hypotheses (complementary cell-state coverage)

*BeatAML2 measured single agents only; no synergy is claimed.*

| pair | pathways | coverage A | coverage B | union | gain |
|---|---|---|---|---|---|
| Cytarabine + Lenvatinib | chemotherapy / RTK | 0.71 | 0.48 | 1.00 | +0.29 |
| Cytarabine + Midostaurin | chemotherapy / FLT3 | 0.71 | 0.45 | 1.00 | +0.29 |
| Cytarabine + Dasatinib | chemotherapy / RTK | 0.71 | 0.81 | 1.00 | +0.19 |
| Cytarabine + Nilotinib | chemotherapy / RTK | 0.71 | 0.84 | 1.00 | +0.16 |
