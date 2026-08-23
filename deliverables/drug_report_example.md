# COMPASS-AML drug prioritisation

*Predicted **ex-vivo** sensitivity from the BeatAML2 functional screen. This is a prioritisation for trial matching or laboratory validation, **not** a treatment recommendation.*

- inhibitors modelled: **118**, reported 59, abstained 59
- genes shared with the model: 14237
- distance from the BeatAML training distribution: 32th percentile
- differentiation axis (primitive - mature): 73th percentile of BeatAML
- cell states scored: 23

## approved in AML  (5 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Quizartinib (AC220) | 0.758 | 0.700 | 0.854 | 1.000 | 1.000 | 0.400 | 0 |
| 2 | Gilteritinib | 0.645 | 0.731 | 0.844 | 0.567 | 0.808 | 0.418 | 1 |
| 3 | Venetoclax | 0.567 | 0.170 | 0.972 | 0.692 | 1.000 | 0.736 | 2 |
| 4 | Midostaurin | 0.397 | 0.760 | 0.798 | 0.171 | 0.000 | 0.476 | 1 |
| 5 | Cytarabine | 0.196 | 0.250 | 0.777 | 0.200 | 0.031 | 0.068 | 1 |

## approved, other indication (off-label in AML)  (20 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Cabozantinib | 0.811 | 0.912 | 0.902 | 0.973 | 0.983 | 0.524 | 1 |
| 2 | Pazopanib (GW786034) | 0.807 | 0.886 | 0.870 | 0.990 | 1.000 | 0.587 | 1 |
| 3 | Regorafenib (BAY 73-4506) | 0.806 | 0.876 | 0.836 | 0.990 | 1.000 | 0.593 | 0 |
| 4 | Lenvatinib | 0.793 | 0.873 | 0.863 | 1.000 | 1.000 | 0.587 | 0 |
| 5 | Dasatinib | 0.785 | 0.895 | 0.915 | 0.990 | 1.000 | 0.556 | 0 |
| 6 | Erlotinib | 0.784 | 0.867 | 0.841 | 0.943 | 1.000 | 0.587 | 0 |
| 7 | Tivozanib (AV-951) | 0.772 | 0.947 | 0.932 | 0.855 | 0.808 | 0.641 | 1 |
| 8 | Sorafenib | 0.768 | 0.864 | 0.919 | 0.839 | 1.000 | 0.502 | 0 |
| 9 | Lapatinib | 0.764 | 0.806 | 0.739 | 0.981 | 0.987 | 0.642 | 1 |
| 10 | Sunitinib | 0.745 | 0.819 | 0.879 | 0.820 | 1.000 | 0.504 | 0 |

## clinical-trial agent  (19 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Linifanib (ABT-869) | 0.762 | 0.910 | 0.914 | 0.990 | 1.000 | 0.511 | 0 |
| 2 | RAF265 (CHIR-265) | 0.755 | 0.866 | 0.840 | 0.986 | 1.000 | 0.721 | 1 |
| 3 | GDC-0941 | 0.743 | 0.781 | 0.826 | 1.000 | 1.000 | 0.677 | 2 |
| 4 | Saracatinib (AZD0530) | 0.741 | 0.863 | 0.864 | 1.000 | 1.000 | 0.574 | 0 |
| 5 | Entospletinib (GS-9973) | 0.730 | 0.961 | 0.885 | 1.000 | 1.000 | 0.356 | 2 |
| 6 | LY-333531 | 0.713 | 0.832 | 0.882 | 0.996 | 1.000 | 0.556 | 1 |
| 7 | Foretinib (XL880) | 0.707 | 0.893 | 0.904 | 0.864 | 0.821 | 0.559 | 1 |
| 8 | Dovitinib (CHIR-258) | 0.688 | 0.871 | 0.921 | 0.686 | 1.000 | 0.514 | 1 |
| 9 | VX-745 | 0.655 | 0.703 | 0.670 | 1.000 | 1.000 | 0.478 | 1 |
| 10 | Pelitinib (EKB-569) | 0.654 | 0.747 | 0.888 | 0.686 | 1.000 | 0.587 | 1 |

## research-only compound  (15 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | NF-kB Activation Inhibitor | 0.696 | 0.938 | 0.862 | 0.712 | 1.000 | 0.809 | 2 |
| 2 | NVP-TAE684 | 0.692 | 0.904 | 0.895 | 0.682 | 0.987 | 0.587 | 1 |
| 3 | PD173955 | 0.687 | 0.819 | 0.800 | 0.990 | 1.000 | 0.534 | 0 |
| 4 | JNJ-28312141 | 0.681 | 0.878 | 0.848 | 0.990 | 1.000 | 0.449 | 2 |
| 5 | SU11274 | 0.681 | 0.792 | 0.791 | 0.849 | 1.000 | 0.575 | 0 |
| 6 | GW-2580 | 0.656 | 0.919 | 0.848 | 0.855 | 0.808 | 0.434 | 2 |
| 7 | Doramapimod (BIRB 796) | 0.645 | 0.785 | 0.852 | 1.000 | 1.000 | 0.444 | 1 |
| 8 | PI-103 | 0.626 | 0.726 | 0.761 | 0.855 | 0.808 | 0.709 | 1 |
| 9 | GSK-1838705A | 0.570 | 0.720 | 0.743 | 0.810 | 0.987 | 0.505 | 2 |
| 10 | PP242 | 0.510 | 0.701 | 0.814 | 0.869 | 0.821 | 0.577 | 0 |

## Combination hypotheses (complementary cell-state coverage)

*BeatAML2 measured single agents only; no synergy is claimed.*

| pair | pathways | coverage A | coverage B | union | gain |
|---|---|---|---|---|---|
| Lenalidomide + Ponatinib (AP24534) | epigenetic / FLT3 | 0.73 | 0.27 | 1.00 | +0.27 |
| Ponatinib (AP24534) + Venetoclax | FLT3 / apoptosis | 0.27 | 0.73 | 1.00 | +0.27 |
