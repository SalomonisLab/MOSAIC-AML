# MOSAIC-Rx drug prioritisation

*Predicted **ex-vivo** sensitivity from the BeatAML2 functional screen. This is a prioritisation for trial matching or laboratory validation, **not** a treatment recommendation.*

- inhibitors modelled: **118**, reported 118, abstained 0
- genes shared with the model: 14237
- distance from the BeatAML training distribution: 32th percentile
- differentiation axis (primitive - mature): 73th percentile of BeatAML
- cell states scored: 23

## approved in AML  (6 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Quizartinib (AC220) | 0.800 | 0.824 | 0.855 | 1.000 | 1.000 | 0.400 | 1 |
| 2 | Gilteritinib | 0.708 | 0.894 | 0.850 | 0.567 | 0.808 | 0.418 | 2 |
| 3 | Venetoclax | 0.597 | 0.317 | 0.966 | 0.692 | 1.000 | 0.736 | 2 |
| 4 | Midostaurin | 0.368 | 0.719 | 0.803 | 0.171 | 0.000 | 0.476 | 1 |
| 5 | Cytarabine | 0.201 | 0.248 | 0.824 | 0.200 | 0.031 | 0.068 | 1 |
| 6 | Azacytidine | 0.156 | 0.345 | 0.646 | 0.000 | 0.000 | 0.091 | 2 |

## approved, other indication (off-label in AML)  (36 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Dasatinib | 0.830 | 0.947 | 0.931 | 0.990 | 1.000 | 0.556 | 0 |
| 2 | Pazopanib (GW786034) | 0.823 | 0.922 | 0.885 | 0.990 | 1.000 | 0.587 | 1 |
| 3 | Cabozantinib | 0.814 | 0.938 | 0.901 | 0.973 | 0.983 | 0.524 | 2 |
| 4 | Lenvatinib | 0.805 | 0.882 | 0.861 | 1.000 | 1.000 | 0.587 | 0 |
| 5 | Regorafenib (BAY 73-4506) | 0.804 | 0.887 | 0.834 | 0.990 | 1.000 | 0.593 | 0 |
| 6 | Tivozanib (AV-951) | 0.786 | 0.976 | 0.941 | 0.855 | 0.808 | 0.641 | 1 |
| 7 | Erlotinib | 0.778 | 0.847 | 0.834 | 0.943 | 1.000 | 0.587 | 0 |
| 8 | Sorafenib | 0.767 | 0.872 | 0.919 | 0.839 | 1.000 | 0.502 | 0 |
| 9 | Lapatinib | 0.757 | 0.794 | 0.734 | 0.981 | 0.987 | 0.642 | 1 |
| 10 | Ibrutinib (PCI-32765) | 0.754 | 0.815 | 0.825 | 1.000 | 1.000 | 0.486 | 0 |

## clinical-trial agent  (36 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | RAF265 (CHIR-265) | 0.781 | 0.885 | 0.831 | 0.986 | 1.000 | 0.721 | 1 |
| 2 | Linifanib (ABT-869) | 0.769 | 0.919 | 0.903 | 0.990 | 1.000 | 0.511 | 1 |
| 3 | Saracatinib (AZD0530) | 0.758 | 0.877 | 0.872 | 1.000 | 1.000 | 0.574 | 0 |
| 4 | Entospletinib (GS-9973) | 0.757 | 0.961 | 0.872 | 1.000 | 1.000 | 0.356 | 2 |
| 5 | GDC-0941 | 0.753 | 0.836 | 0.825 | 1.000 | 1.000 | 0.677 | 2 |
| 6 | LY-333531 | 0.741 | 0.838 | 0.871 | 0.996 | 1.000 | 0.556 | 1 |
| 7 | Foretinib (XL880) | 0.715 | 0.930 | 0.908 | 0.864 | 0.821 | 0.559 | 1 |
| 8 | Dovitinib (CHIR-258) | 0.708 | 0.889 | 0.917 | 0.686 | 1.000 | 0.514 | 1 |
| 9 | Pelitinib (EKB-569) | 0.699 | 0.845 | 0.893 | 0.686 | 1.000 | 0.587 | 1 |
| 10 | SNS-032 (BMS-387032) | 0.674 | 0.849 | 0.801 | 0.990 | 1.000 | 0.192 | 1 |

## research-only compound  (40 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | NF-kB Activation Inhibitor | 0.722 | 0.932 | 0.865 | 0.712 | 1.000 | 0.809 | 2 |
| 2 | NVP-TAE684 | 0.717 | 0.931 | 0.907 | 0.682 | 0.987 | 0.587 | 1 |
| 3 | PD173955 | 0.716 | 0.855 | 0.794 | 0.990 | 1.000 | 0.534 | 0 |
| 4 | JNJ-28312141 | 0.693 | 0.901 | 0.842 | 0.990 | 1.000 | 0.449 | 2 |
| 5 | ABT-737 | 0.669 | 0.789 | 0.763 | 0.701 | 1.000 | 0.782 | 1 |
| 6 | GW-2580 | 0.669 | 0.937 | 0.863 | 0.855 | 0.808 | 0.434 | 2 |
| 7 | SU11274 | 0.658 | 0.753 | 0.809 | 0.849 | 1.000 | 0.575 | 0 |
| 8 | PI-103 | 0.640 | 0.795 | 0.757 | 0.855 | 0.808 | 0.709 | 1 |
| 9 | Doramapimod (BIRB 796) | 0.625 | 0.740 | 0.841 | 1.000 | 1.000 | 0.444 | 1 |
| 10 | KU-55933 | 0.594 | 0.668 | 0.737 | 0.986 | 1.000 | 0.496 | 2 |

## Combination hypotheses (complementary cell-state coverage)

*BeatAML2 measured single agents only; no synergy is claimed.*

| pair | pathways | coverage A | coverage B | union | gain |
|---|---|---|---|---|---|
| Azacytidine + Neratinib (HKI-272) | epigenetic / RTK | 0.27 | 0.61 | 0.88 | +0.27 |
| Azacytidine + Palbociclib | epigenetic / cell_cycle | 0.27 | 0.73 | 1.00 | +0.27 |
| Bortezomib (Velcade) + Lenalidomide | proteostasis / epigenetic | 0.27 | 0.73 | 1.00 | +0.27 |
| Bortezomib (Velcade) + Neratinib (HKI-272) | proteostasis / RTK | 0.27 | 0.61 | 0.88 | +0.27 |
| Bortezomib (Velcade) + Palbociclib | proteostasis / cell_cycle | 0.27 | 0.73 | 1.00 | +0.27 |
| Idelalisib + Lenalidomide | PI3K_AKT_mTOR / epigenetic | 0.27 | 0.73 | 1.00 | +0.27 |
| Idelalisib + Neratinib (HKI-272) | PI3K_AKT_mTOR / RTK | 0.27 | 0.61 | 0.88 | +0.27 |
| Idelalisib + Palbociclib | PI3K_AKT_mTOR / cell_cycle | 0.27 | 0.73 | 1.00 | +0.27 |
