# MOSAIC-Rx drug prioritisation

*Predicted **ex-vivo** sensitivity from the BeatAML2 functional screen. This is a prioritisation for trial matching or laboratory validation, **not** a treatment recommendation.*

- inhibitors modelled: **118**, reported 118, abstained 0
- genes shared with the model: 14237
- distance from the BeatAML training distribution: 69th percentile
- differentiation axis (primitive - mature): 54th percentile of BeatAML
- cell states scored: 30

## approved in AML  (6 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Quizartinib (AC220) | 0.792 | 0.902 | 0.855 | 0.995 | 1.000 | 0.651 | 0 |
| 2 | Gilteritinib | 0.756 | 0.919 | 0.850 | 0.652 | 1.000 | 0.747 | 0 |
| 3 | Midostaurin | 0.741 | 0.844 | 0.803 | 0.879 | 1.000 | 0.632 | 0 |
| 4 | Cytarabine | 0.702 | 0.831 | 0.824 | 0.936 | 1.000 | 0.521 | 1 |
| 5 | Venetoclax | 0.411 | 0.021 | 0.966 | 0.495 | 1.000 | 0.536 | 2 |
| 6 | Azacytidine | 0.222 | 0.321 | 0.646 | 0.056 | 0.107 | 0.497 | 1 |

## approved, other indication (off-label in AML)  (36 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Tivozanib (AV-951) | 0.779 | 0.986 | 0.941 | 0.991 | 1.000 | 0.592 | 0 |
| 2 | Dasatinib | 0.777 | 0.965 | 0.931 | 0.989 | 1.000 | 0.646 | 0 |
| 3 | Ponatinib (AP24534) | 0.770 | 0.981 | 0.952 | 0.907 | 1.000 | 0.643 | 0 |
| 4 | Bosutinib (SKI-606) | 0.768 | 0.926 | 0.862 | 0.991 | 1.000 | 0.725 | 0 |
| 5 | Sunitinib | 0.764 | 0.946 | 0.891 | 1.000 | 1.000 | 0.642 | 0 |
| 6 | Pazopanib (GW786034) | 0.760 | 0.948 | 0.885 | 1.000 | 1.000 | 0.600 | 0 |
| 7 | Sorafenib | 0.750 | 0.964 | 0.919 | 0.894 | 1.000 | 0.606 | 0 |
| 8 | Ibrutinib (PCI-32765) | 0.743 | 0.911 | 0.825 | 1.000 | 1.000 | 0.643 | 0 |
| 9 | Gefitinib | 0.741 | 0.878 | 0.794 | 0.967 | 1.000 | 0.770 | 0 |
| 10 | Regorafenib (BAY 73-4506) | 0.738 | 0.913 | 0.834 | 0.989 | 1.000 | 0.610 | 0 |

## clinical-trial agent  (36 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Saracatinib (AZD0530) | 0.742 | 0.940 | 0.872 | 1.000 | 1.000 | 0.769 | 0 |
| 2 | Linifanib (ABT-869) | 0.721 | 0.948 | 0.903 | 0.980 | 1.000 | 0.636 | 0 |
| 3 | KW-2449 | 0.717 | 0.957 | 0.913 | 0.894 | 1.000 | 0.689 | 0 |
| 4 | Entospletinib (GS-9973) | 0.716 | 0.952 | 0.872 | 0.864 | 1.000 | 0.752 | 0 |
| 5 | Crenolanib | 0.659 | 0.902 | 0.879 | 0.736 | 1.000 | 0.639 | 0 |
| 6 | GDC-0941 | 0.658 | 0.806 | 0.825 | 0.971 | 0.985 | 0.706 | 1 |
| 7 | LY-333531 | 0.655 | 0.907 | 0.871 | 0.655 | 1.000 | 0.707 | 0 |
| 8 | Dovitinib (CHIR-258) | 0.653 | 0.915 | 0.917 | 0.634 | 1.000 | 0.658 | 0 |
| 9 | Foretinib (XL880) | 0.640 | 0.961 | 0.908 | 0.988 | 1.000 | 0.637 | 0 |
| 10 | Pelitinib (EKB-569) | 0.632 | 0.854 | 0.893 | 0.585 | 1.000 | 0.770 | 0 |

## research-only compound  (40 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | NVP-TAE684 | 0.725 | 0.968 | 0.907 | 0.919 | 1.000 | 0.770 | 0 |
| 2 | PRT062607 | 0.685 | 0.877 | 0.804 | 0.995 | 1.000 | 0.752 | 0 |
| 3 | PD173955 | 0.680 | 0.874 | 0.794 | 0.991 | 1.000 | 0.725 | 0 |
| 4 | JNJ-28312141 | 0.643 | 0.872 | 0.842 | 0.941 | 1.000 | 0.740 | 1 |
| 5 | Doramapimod (BIRB 796) | 0.634 | 0.909 | 0.841 | 1.000 | 1.000 | 0.536 | 1 |
| 6 | NVP-ADW742 | 0.602 | 0.805 | 0.723 | 0.974 | 1.000 | 0.707 | 1 |
| 7 | JAK Inhibitor I | 0.591 | 0.711 | 0.700 | 1.000 | 1.000 | 0.642 | 0 |
| 8 | SU11274 | 0.586 | 0.871 | 0.809 | 1.000 | 1.000 | 0.712 | 0 |
| 9 | TG100-115 | 0.550 | 0.703 | 0.689 | 0.930 | 0.878 | 0.715 | 0 |
| 10 | PP242 | 0.529 | 0.844 | 0.816 | 0.939 | 0.884 | 0.660 | 0 |

## Combination hypotheses (complementary cell-state coverage)

*BeatAML2 measured single agents only; no synergy is claimed.*

| pair | pathways | coverage A | coverage B | union | gain |
|---|---|---|---|---|---|
| Neratinib (HKI-272) + Panobinostat | RTK / epigenetic | 0.49 | 0.51 | 1.00 | +0.49 |
| Erlotinib + Panobinostat | RTK / epigenetic | 0.64 | 0.51 | 1.00 | +0.36 |
| Palbociclib + Panobinostat | cell_cycle / epigenetic | 0.64 | 0.51 | 1.00 | +0.36 |
| Afatinib (BIBW-2992) + Panobinostat | RTK / epigenetic | 0.84 | 0.51 | 1.00 | +0.15 |
