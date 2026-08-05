# MOSAIC-Rx drug prioritisation

*Predicted **ex-vivo** sensitivity from the BeatAML2 functional screen. This is a prioritisation for trial matching or laboratory validation, **not** a treatment recommendation.*

- inhibitors modelled: **118**, reported 117, abstained 1
- genes shared with the model: 14237
- distance from the BeatAML training distribution: 97th percentile
- differentiation axis (primitive - mature): 9th percentile of BeatAML
- cell states scored: 15

## approved in AML  (6 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Midostaurin | 0.468 | 0.764 | 0.803 | 0.763 | 0.792 | 0.301 | 3 |
| 2 | Quizartinib (AC220) | 0.440 | 0.565 | 0.855 | 0.745 | 1.000 | 0.292 | 1 |
| 3 | Gilteritinib | 0.417 | 0.432 | 0.850 | 0.753 | 1.000 | 0.424 | 2 |
| 4 | Cytarabine | 0.000 | 0.327 | 0.824 | 0.040 | 0.208 | 0.000 | 2 |
| 5 | Venetoclax | -0.072 | 0.003 | 0.966 | 0.000 | 0.000 | 0.000 | 3 |
| 6 | Azacytidine | -0.100 | 0.253 | 0.646 | 0.000 | 0.000 | 0.000 | 3 |

## approved, other indication (off-label in AML)  (36 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Tivozanib (AV-951) | 0.671 | 0.984 | 0.941 | 1.000 | 1.000 | 0.612 | 2 |
| 2 | Trametinib (GSK1120212) | 0.652 | 0.962 | 0.903 | 1.000 | 1.000 | 0.617 | 1 |
| 3 | Dasatinib | 0.640 | 0.968 | 0.931 | 1.000 | 1.000 | 0.477 | 2 |
| 4 | Pazopanib (GW786034) | 0.628 | 0.937 | 0.885 | 1.000 | 1.000 | 0.518 | 3 |
| 5 | Regorafenib (BAY 73-4506) | 0.627 | 0.923 | 0.834 | 1.000 | 1.000 | 0.582 | 2 |
| 6 | Selumetinib (AZD6244) | 0.592 | 0.952 | 0.886 | 0.994 | 0.792 | 0.617 | 1 |
| 7 | Imatinib | 0.538 | 0.773 | 0.707 | 1.000 | 1.000 | 0.491 | 2 |
| 8 | CYT387 | 0.520 | 0.755 | 0.703 | 0.995 | 1.000 | 0.452 | 1 |
| 9 | Bosutinib (SKI-606) | 0.517 | 0.686 | 0.862 | 0.748 | 1.000 | 0.742 | 2 |
| 10 | Lapatinib | 0.494 | 0.720 | 0.734 | 0.970 | 1.000 | 0.377 | 2 |

## clinical-trial agent  (36 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Saracatinib (AZD0530) | 0.617 | 0.925 | 0.872 | 1.000 | 1.000 | 0.741 | 1 |
| 2 | CI-1040 (PD184352) | 0.578 | 0.897 | 0.819 | 1.000 | 1.000 | 0.617 | 1 |
| 3 | GDC-0941 | 0.535 | 0.824 | 0.825 | 0.995 | 1.000 | 0.540 | 4 |
| 4 | Motesanib (AMG-706) | 0.493 | 0.862 | 0.824 | 0.999 | 0.792 | 0.527 | 1 |
| 5 | Flavopiridol | 0.374 | 0.737 | 0.695 | 0.999 | 0.792 | 0.185 | 3 |
| 6 | 17-AAG (Tanespimycin) | 0.373 | 0.958 | 0.894 | 0.992 | 0.415 | 0.103 | 2 |
| 7 | SNS-032 (BMS-387032) | 0.359 | 0.647 | 0.801 | 0.949 | 0.792 | 0.298 | 3 |
| 8 | Roscovitine (CYC-202) | 0.358 | 0.580 | 0.603 | 0.986 | 1.000 | 0.241 | 3 |
| 9 | Foretinib (XL880) | 0.333 | 0.726 | 0.908 | 0.750 | 0.415 | 0.660 | 2 |
| 10 | MK-2206 | 0.325 | 0.784 | 0.767 | 0.994 | 0.792 | 0.297 | 3 |

## research-only compound  (39 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | PRT062607 | 0.554 | 0.861 | 0.804 | 0.995 | 1.000 | 0.676 | 1 |
| 2 | PD173955 | 0.553 | 0.825 | 0.794 | 1.000 | 1.000 | 0.742 | 2 |
| 3 | GW-2580 | 0.531 | 0.941 | 0.863 | 0.983 | 0.792 | 0.630 | 2 |
| 4 | PI-103 | 0.501 | 0.839 | 0.757 | 1.000 | 1.000 | 0.448 | 3 |
| 5 | JNJ-28312141 | 0.496 | 0.869 | 0.842 | 0.984 | 1.000 | 0.493 | 3 |
| 6 | Doramapimod (BIRB 796) | 0.483 | 0.898 | 0.841 | 0.995 | 1.000 | 0.330 | 3 |
| 7 | NVP-ADW742 | 0.441 | 0.733 | 0.723 | 0.986 | 1.000 | 0.565 | 2 |
| 8 | PP242 | 0.410 | 0.890 | 0.816 | 0.995 | 1.000 | 0.265 | 3 |
| 9 | TG100-115 | 0.400 | 0.710 | 0.689 | 0.994 | 0.792 | 0.543 | 2 |
| 10 | JNJ-7706621 | 0.354 | 0.649 | 0.663 | 0.991 | 1.000 | 0.096 | 3 |
