# MOSAIC-Rx drug prioritisation

*Predicted **ex-vivo** sensitivity from the BeatAML2 functional screen. This is a prioritisation for trial matching or laboratory validation, **not** a treatment recommendation.*

- inhibitors modelled: **118**, reported 118, abstained 0
- genes shared with the model: 14237
- distance from the BeatAML training distribution: 81th percentile
- differentiation axis (primitive - mature): 50th percentile of BeatAML
- cell states scored: 46

## approved in AML  (6 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Venetoclax | 0.651 | 0.950 | 0.966 | 0.831 | 1.000 | 0.440 | 2 |
| 2 | Midostaurin | 0.415 | 0.562 | 0.803 | 0.679 | 0.635 | 0.429 | 0 |
| 3 | Gilteritinib | 0.253 | 0.270 | 0.850 | 0.393 | 0.309 | 0.353 | 0 |
| 4 | Cytarabine | 0.206 | 0.561 | 0.824 | 0.266 | 0.000 | 0.565 | 2 |
| 5 | Azacytidine | 0.196 | 0.447 | 0.646 | 0.048 | 0.000 | 0.523 | 2 |
| 6 | Quizartinib (AC220) | 0.126 | 0.143 | 0.855 | 0.019 | 0.000 | 0.437 | 0 |

## approved, other indication (off-label in AML)  (36 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Ponatinib (AP24534) | 0.611 | 0.863 | 0.952 | 0.814 | 0.944 | 0.449 | 1 |
| 2 | Bosutinib (SKI-606) | 0.609 | 0.869 | 0.862 | 0.910 | 1.000 | 0.297 | 1 |
| 3 | Neratinib (HKI-272) | 0.550 | 0.765 | 0.716 | 0.913 | 1.000 | 0.295 | 2 |
| 4 | Palbociclib | 0.509 | 0.683 | 0.782 | 0.850 | 0.963 | 0.604 | 0 |
| 5 | Lenalidomide | 0.505 | 0.607 | 0.726 | 0.819 | 1.000 | 0.518 | 0 |
| 6 | Bortezomib (Velcade) | 0.486 | 0.699 | 0.667 | 0.997 | 1.000 | 0.553 | 1 |
| 7 | Axitinib (AG-013736) | 0.469 | 0.639 | 0.862 | 0.884 | 0.819 | 0.331 | 0 |
| 8 | Trametinib (GSK1120212) | 0.389 | 0.528 | 0.903 | 0.793 | 0.880 | 0.125 | 1 |
| 9 | Selumetinib (AZD6244) | 0.374 | 0.488 | 0.886 | 0.769 | 0.847 | 0.125 | 0 |
| 10 | Selinexor | 0.361 | 0.501 | 0.550 | 0.610 | 0.865 | 0.456 | 2 |

## clinical-trial agent  (36 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Volasertib (BI-6727) | 0.559 | 0.690 | 0.668 | 0.931 | 1.000 | 0.804 | 0 |
| 2 | OTX-015 | 0.543 | 0.706 | 0.779 | 0.916 | 1.000 | 0.627 | 0 |
| 3 | Cediranib (AZD2171) | 0.489 | 0.676 | 0.748 | 0.988 | 0.994 | 0.331 | 0 |
| 4 | Motesanib (AMG-706) | 0.445 | 0.600 | 0.824 | 0.961 | 0.939 | 0.309 | 0 |
| 5 | Elesclomol | 0.439 | 0.785 | 0.800 | 0.951 | 0.994 | 0.024 | 2 |
| 6 | GSK690693 | 0.431 | 0.556 | 0.658 | 0.910 | 1.000 | 0.402 | 0 |
| 7 | MK-2206 | 0.408 | 0.462 | 0.767 | 0.917 | 0.994 | 0.402 | 0 |
| 8 | AZD1480 | 0.378 | 0.580 | 0.680 | 0.843 | 0.870 | 0.281 | 0 |
| 9 | Masitinib (AB-1010) | 0.288 | 0.653 | 0.684 | 0.718 | 0.467 | 0.285 | 0 |
| 10 | Dovitinib (CHIR-258) | 0.278 | 0.519 | 0.917 | 0.629 | 0.522 | 0.412 | 0 |

## research-only compound  (40 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Nutlin 3a | 0.499 | 0.765 | 0.786 | 0.856 | 0.898 | 0.518 | 1 |
| 2 | JQ1 | 0.475 | 0.724 | 0.763 | 0.824 | 0.831 | 0.627 | 1 |
| 3 | ABT-737 | 0.418 | 0.589 | 0.763 | 0.724 | 0.870 | 0.669 | 1 |
| 4 | DBZ | 0.381 | 0.589 | 0.667 | 0.831 | 1.000 | 0.380 | 1 |
| 5 | JAK Inhibitor I | 0.362 | 0.390 | 0.700 | 0.870 | 0.939 | 0.446 | 0 |
| 6 | BMS-345541 | 0.330 | 0.554 | 0.556 | 0.903 | 0.957 | 0.223 | 2 |
| 7 | Bay 11-7085 | 0.317 | 0.549 | 0.705 | 0.827 | 0.963 | 0.148 | 3 |
| 8 | GW-2580 | 0.221 | 0.441 | 0.863 | 0.477 | 0.198 | 0.564 | 1 |
| 9 | Artemisinin | 0.219 | 0.573 | 0.665 | 0.572 | 0.388 | — | 1 |
| 10 | PRT062607 | 0.199 | 0.418 | 0.804 | 0.573 | 0.326 | 0.294 | 0 |

## Combination hypotheses (complementary cell-state coverage)

*BeatAML2 measured single agents only; no synergy is claimed.*

| pair | pathways | coverage A | coverage B | union | gain |
|---|---|---|---|---|---|
| CYT387 + Cytarabine | JAK_STAT / chemotherapy | 0.27 | 0.62 | 0.85 | +0.23 |
| CYT387 + Entrectinib | JAK_STAT / RTK | 0.27 | 0.70 | 0.93 | +0.23 |
| Cytarabine + Gilteritinib | chemotherapy / FLT3 | 0.62 | 0.51 | 0.84 | +0.22 |
| Cytarabine + Dasatinib | chemotherapy / RTK | 0.62 | 0.43 | 0.83 | +0.21 |
| Cytarabine + Lapatinib | chemotherapy / RTK | 0.62 | 0.51 | 0.82 | +0.20 |
| Afatinib (BIBW-2992) + Cytarabine | RTK / chemotherapy | 0.29 | 0.62 | 0.81 | +0.19 |
| Cabozantinib + Cytarabine | RTK / chemotherapy | 0.59 | 0.62 | 0.81 | +0.19 |
| Cytarabine + Gefitinib | chemotherapy / RTK | 0.62 | 0.19 | 0.81 | +0.19 |
