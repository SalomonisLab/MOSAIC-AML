# MOSAIC-Rx drug prioritisation

*Predicted **ex-vivo** sensitivity from the BeatAML2 functional screen. This is a prioritisation for trial matching or laboratory validation, **not** a treatment recommendation.*

- inhibitors modelled: **118**, reported 118, abstained 0
- genes shared with the model: 14237
- distance from the BeatAML training distribution: 44th percentile
- differentiation axis (primitive - mature): 55th percentile of BeatAML
- cell states scored: 38

## approved in AML  (6 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Venetoclax | 0.665 | 0.506 | 0.966 | 0.792 | 1.000 | 0.562 | 1 |
| 2 | Azacytidine | 0.655 | 0.563 | 0.646 | 0.899 | 1.000 | 0.451 | 2 |
| 3 | Quizartinib (AC220) | 0.607 | 0.565 | 0.855 | 0.654 | 0.787 | 0.584 | 0 |
| 4 | Cytarabine | 0.554 | 0.336 | 0.824 | 0.749 | 0.916 | 0.360 | 1 |
| 5 | Midostaurin | 0.386 | 0.361 | 0.803 | 0.215 | 0.237 | 0.551 | 1 |
| 6 | Gilteritinib | 0.279 | 0.206 | 0.850 | 0.060 | 0.052 | 0.460 | 0 |

## approved, other indication (off-label in AML)  (36 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Ponatinib (AP24534) | 0.790 | 0.979 | 0.952 | 0.869 | 0.921 | 0.473 | 2 |
| 2 | Sorafenib | 0.698 | 0.801 | 0.919 | 0.774 | 0.916 | 0.491 | 0 |
| 3 | Tivozanib (AV-951) | 0.688 | 0.748 | 0.941 | 0.864 | 0.883 | 0.502 | 0 |
| 4 | Vargetef | 0.681 | 0.724 | 0.827 | 0.862 | 0.943 | 0.505 | 0 |
| 5 | Lapatinib | 0.672 | 0.725 | 0.734 | 0.859 | 0.994 | 0.430 | 1 |
| 6 | Crizotinib (PF-2341066) | 0.668 | 0.782 | 0.802 | 0.929 | 0.934 | 0.242 | 1 |
| 7 | Lenvatinib | 0.656 | 0.713 | 0.861 | 0.862 | 0.847 | 0.515 | 0 |
| 8 | Axitinib (AG-013736) | 0.650 | 0.719 | 0.862 | 0.833 | 0.798 | 0.549 | 0 |
| 9 | Afatinib (BIBW-2992) | 0.640 | 0.688 | 0.718 | 0.847 | 0.943 | 0.430 | 0 |
| 10 | Neratinib (HKI-272) | 0.634 | 0.668 | 0.716 | 0.836 | 0.962 | 0.430 | 1 |

## clinical-trial agent  (36 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Pelitinib (EKB-569) | 0.693 | 0.858 | 0.893 | 0.848 | 1.000 | 0.367 | 2 |
| 2 | Volasertib (BI-6727) | 0.685 | 0.724 | 0.668 | 0.943 | 0.949 | 0.796 | 0 |
| 3 | OTX-015 | 0.684 | 0.836 | 0.779 | 0.989 | 1.000 | 0.319 | 1 |
| 4 | KW-2449 | 0.666 | 0.806 | 0.913 | 0.841 | 0.934 | 0.444 | 0 |
| 5 | Alisertib (MLN8237) | 0.647 | 0.733 | 0.728 | 0.819 | 0.943 | 0.632 | 1 |
| 6 | Masitinib (AB-1010) | 0.640 | 0.701 | 0.684 | 0.951 | 1.000 | 0.470 | 0 |
| 7 | AZD1152-HQPA (AZD2811) | 0.628 | 0.655 | 0.824 | 0.744 | 0.923 | 0.761 | 1 |
| 8 | Canertinib (CI-1033) | 0.615 | 0.677 | 0.773 | 0.848 | 1.000 | 0.430 | 1 |
| 9 | MGCD-265 | 0.570 | 0.736 | 0.843 | 0.742 | 0.839 | 0.337 | 1 |
| 10 | Roscovitine (CYC-202) | 0.566 | 0.583 | 0.603 | 0.822 | 0.967 | 0.592 | 1 |

## research-only compound  (40 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | GW-2580 | 0.748 | 0.903 | 0.863 | 0.975 | 1.000 | 0.596 | 0 |
| 2 | JQ1 | 0.661 | 0.827 | 0.763 | 0.979 | 1.000 | 0.319 | 1 |
| 3 | AGI-6780 | 0.639 | 0.753 | 0.773 | 0.760 | 1.000 | 0.598 | 1 |
| 4 | PHA-665752 | 0.623 | 0.828 | 0.798 | 0.808 | 1.000 | 0.242 | 1 |
| 5 | KI20227 | 0.601 | 0.776 | 0.808 | 0.831 | 0.883 | 0.589 | 2 |
| 6 | KU-55933 | 0.597 | 0.760 | 0.737 | 0.887 | 1.000 | 0.406 | 1 |
| 7 | Artemisinin | 0.580 | 0.633 | 0.665 | 0.836 | 0.934 | — | 2 |
| 8 | PI-103 | 0.569 | 0.775 | 0.757 | 0.889 | 0.878 | 0.179 | 1 |
| 9 | Bay 11-7085 | 0.568 | 0.682 | 0.705 | 0.917 | 1.000 | 0.422 | 1 |
| 10 | PD173955 | 0.567 | 0.702 | 0.794 | 0.808 | 0.943 | 0.294 | 1 |

## Combination hypotheses (complementary cell-state coverage)

*BeatAML2 measured single agents only; no synergy is claimed.*

| pair | pathways | coverage A | coverage B | union | gain |
|---|---|---|---|---|---|
| Gefitinib + Midostaurin | RTK / FLT3 | 0.21 | 0.20 | 0.37 | +0.16 |
