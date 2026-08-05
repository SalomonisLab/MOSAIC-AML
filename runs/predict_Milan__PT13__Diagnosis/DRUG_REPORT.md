# MOSAIC-Rx drug prioritisation

*Predicted **ex-vivo** sensitivity from the BeatAML2 functional screen. This is a prioritisation for trial matching or laboratory validation, **not** a treatment recommendation.*

- inhibitors modelled: **118**, reported 118, abstained 0
- genes shared with the model: 14237
- distance from the BeatAML training distribution: 8th percentile
- differentiation axis (primitive - mature): 95th percentile of BeatAML
- cell states scored: 17

## approved in AML  (6 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Venetoclax | 0.941 | 0.977 | 0.966 | 0.972 | 1.000 | 0.854 | 1 |
| 2 | Quizartinib (AC220) | 0.830 | 0.868 | 0.855 | 0.977 | 1.000 | 0.514 | 0 |
| 3 | Cytarabine | 0.790 | 0.875 | 0.824 | 0.887 | 1.000 | 0.381 | 2 |
| 4 | Gilteritinib | 0.517 | 0.440 | 0.850 | 0.551 | 0.592 | 0.487 | 0 |
| 5 | Midostaurin | 0.436 | 0.332 | 0.803 | 0.421 | 0.463 | 0.450 | 0 |
| 6 | Azacytidine | 0.280 | 0.380 | 0.646 | 0.014 | 0.016 | 0.393 | 1 |

## approved, other indication (off-label in AML)  (36 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Lenvatinib | 0.812 | 0.911 | 0.861 | 0.987 | 1.000 | 0.590 | 0 |
| 2 | Ponatinib (AP24534) | 0.810 | 0.958 | 0.952 | 0.901 | 1.000 | 0.494 | 1 |
| 3 | Sorafenib | 0.809 | 0.956 | 0.919 | 0.900 | 1.000 | 0.512 | 1 |
| 4 | Cabozantinib | 0.801 | 0.929 | 0.901 | 0.893 | 1.000 | 0.536 | 0 |
| 5 | Sunitinib | 0.796 | 0.914 | 0.891 | 0.893 | 1.000 | 0.551 | 0 |
| 6 | Erlotinib | 0.795 | 0.861 | 0.834 | 0.957 | 1.000 | 0.654 | 0 |
| 7 | Pazopanib (GW786034) | 0.774 | 0.839 | 0.885 | 0.901 | 1.000 | 0.601 | 0 |
| 8 | Gefitinib | 0.773 | 0.802 | 0.794 | 0.984 | 1.000 | 0.654 | 0 |
| 9 | Regorafenib (BAY 73-4506) | 0.768 | 0.853 | 0.834 | 0.900 | 1.000 | 0.560 | 0 |
| 10 | Vargetef | 0.761 | 0.837 | 0.827 | 0.878 | 1.000 | 0.594 | 0 |

## clinical-trial agent  (36 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Dovitinib (CHIR-258) | 0.797 | 0.935 | 0.917 | 0.969 | 1.000 | 0.665 | 0 |
| 2 | Pelitinib (EKB-569) | 0.783 | 0.906 | 0.893 | 0.992 | 1.000 | 0.654 | 0 |
| 3 | LY-333531 | 0.755 | 0.886 | 0.871 | 0.979 | 1.000 | 0.566 | 0 |
| 4 | KW-2449 | 0.750 | 0.912 | 0.913 | 0.893 | 1.000 | 0.526 | 0 |
| 5 | Entospletinib (GS-9973) | 0.742 | 0.892 | 0.872 | 0.980 | 1.000 | 0.469 | 1 |
| 6 | Foretinib (XL880) | 0.742 | 0.881 | 0.908 | 0.886 | 1.000 | 0.544 | 0 |
| 7 | Linifanib (ABT-869) | 0.704 | 0.798 | 0.903 | 0.887 | 1.000 | 0.539 | 0 |
| 8 | Crenolanib | 0.700 | 0.816 | 0.879 | 0.890 | 1.000 | 0.484 | 0 |
| 9 | Saracatinib (AZD0530) | 0.696 | 0.745 | 0.872 | 0.901 | 1.000 | 0.630 | 0 |
| 10 | SNS-032 (BMS-387032) | 0.652 | 0.763 | 0.801 | 0.901 | 1.000 | 0.365 | 1 |

## research-only compound  (40 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | NVP-TAE684 | 0.757 | 0.876 | 0.907 | 0.992 | 1.000 | 0.654 | 0 |
| 2 | Nutlin 3a | 0.739 | 0.849 | 0.786 | 0.992 | 1.000 | 0.721 | 2 |
| 3 | NF-kB Activation Inhibitor | 0.713 | 0.932 | 0.865 | 0.970 | 0.984 | 0.525 | 2 |
| 4 | ABT-737 | 0.701 | 0.752 | 0.763 | 0.951 | 0.984 | 0.812 | 1 |
| 5 | SU11274 | 0.691 | 0.824 | 0.809 | 0.980 | 1.000 | 0.475 | 0 |
| 6 | PHA-665752 | 0.663 | 0.794 | 0.798 | 0.893 | 1.000 | 0.475 | 0 |
| 7 | PD173955 | 0.641 | 0.676 | 0.794 | 0.901 | 1.000 | 0.618 | 0 |
| 8 | JNJ-28312141 | 0.640 | 0.756 | 0.842 | 0.894 | 1.000 | 0.588 | 1 |
| 9 | A-674563 | 0.631 | 0.694 | 0.777 | 0.878 | 1.000 | 0.560 | 1 |
| 10 | PRT062607 | 0.617 | 0.651 | 0.804 | 0.973 | 1.000 | 0.469 | 0 |
