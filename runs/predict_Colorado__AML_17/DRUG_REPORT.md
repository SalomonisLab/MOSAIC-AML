# MOSAIC-Rx drug prioritisation

*Predicted **ex-vivo** sensitivity from the BeatAML2 functional screen. This is a prioritisation for trial matching or laboratory validation, **not** a treatment recommendation.*

- inhibitors modelled: **118**, reported 118, abstained 0
- genes shared with the model: 14237
- distance from the BeatAML training distribution: 46th percentile
- differentiation axis (primitive - mature): 9th percentile of BeatAML
- cell states scored: 22

## approved in AML  (6 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Midostaurin | 0.735 | 0.715 | 0.803 | 0.859 | 1.000 | 0.477 | 0 |
| 2 | Quizartinib (AC220) | 0.693 | 0.565 | 0.855 | 0.977 | 1.000 | 0.449 | 0 |
| 3 | Gilteritinib | 0.477 | 0.335 | 0.850 | 0.121 | 1.000 | 0.435 | 0 |
| 4 | Venetoclax | 0.259 | 0.002 | 0.966 | 0.049 | 1.000 | 0.000 | 2 |
| 5 | Azacytidine | 0.196 | 0.236 | 0.646 | 0.034 | 0.107 | 0.198 | 2 |
| 6 | Cytarabine | 0.187 | 0.146 | 0.824 | 0.090 | 0.107 | 0.223 | 1 |

## approved, other indication (off-label in AML)  (36 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Rapamycin | 0.831 | 0.947 | 0.922 | 0.997 | 0.893 | 0.727 | 1 |
| 2 | Tivozanib (AV-951) | 0.824 | 0.956 | 0.941 | 0.998 | 1.000 | 0.483 | 1 |
| 3 | Ibrutinib (PCI-32765) | 0.815 | 0.882 | 0.825 | 1.000 | 1.000 | 0.689 | 0 |
| 4 | Dasatinib | 0.804 | 0.975 | 0.931 | 0.994 | 0.893 | 0.495 | 1 |
| 5 | Lenvatinib | 0.790 | 0.915 | 0.861 | 0.992 | 1.000 | 0.440 | 1 |
| 6 | Pazopanib (GW786034) | 0.773 | 0.861 | 0.885 | 0.998 | 1.000 | 0.441 | 1 |
| 7 | CYT387 | 0.763 | 0.783 | 0.703 | 1.000 | 1.000 | 0.681 | 0 |
| 8 | Entrectinib | 0.760 | 0.874 | 0.771 | 0.998 | 1.000 | 0.435 | 2 |
| 9 | Ruxolitinib (INCB018424) | 0.756 | 0.785 | 0.713 | 0.998 | 1.000 | 0.648 | 0 |
| 10 | Cabozantinib | 0.755 | 0.812 | 0.901 | 0.990 | 1.000 | 0.455 | 0 |

## clinical-trial agent  (36 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | 17-AAG (Tanespimycin) | 0.813 | 0.940 | 0.894 | 0.998 | 1.000 | 0.734 | 0 |
| 2 | Saracatinib (AZD0530) | 0.799 | 0.916 | 0.872 | 1.000 | 1.000 | 0.706 | 0 |
| 3 | Entospletinib (GS-9973) | 0.788 | 0.860 | 0.872 | 0.994 | 1.000 | 0.838 | 0 |
| 4 | INK-128 | 0.758 | 0.836 | 0.793 | 1.000 | 1.000 | 0.727 | 1 |
| 5 | GDC-0941 | 0.751 | 0.898 | 0.825 | 0.997 | 0.893 | 0.676 | 2 |
| 6 | Foretinib (XL880) | 0.748 | 0.902 | 0.908 | 0.998 | 1.000 | 0.418 | 1 |
| 7 | MGCD-265 | 0.717 | 0.846 | 0.843 | 0.994 | 1.000 | 0.428 | 2 |
| 8 | CI-1040 (PD184352) | 0.712 | 0.881 | 0.819 | 0.997 | 0.893 | 0.490 | 0 |
| 9 | AZD1480 | 0.701 | 0.717 | 0.680 | 1.000 | 1.000 | 0.767 | 0 |
| 10 | OTX-015 | 0.695 | 0.862 | 0.779 | 0.997 | 0.893 | 0.467 | 0 |

## research-only compound  (40 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | PRT062607 | 0.778 | 0.882 | 0.804 | 1.000 | 1.000 | 0.838 | 0 |
| 2 | PP242 | 0.734 | 0.878 | 0.816 | 0.994 | 0.893 | 0.727 | 1 |
| 3 | Doramapimod (BIRB 796) | 0.731 | 0.913 | 0.841 | 1.000 | 1.000 | 0.643 | 1 |
| 4 | PD173955 | 0.717 | 0.785 | 0.794 | 0.998 | 1.000 | 0.704 | 0 |
| 5 | JAK Inhibitor I | 0.680 | 0.747 | 0.700 | 1.000 | 1.000 | 0.624 | 0 |
| 6 | TG100-115 | 0.673 | 0.683 | 0.689 | 0.998 | 1.000 | 0.761 | 1 |
| 7 | JNJ-28312141 | 0.666 | 0.773 | 0.842 | 0.963 | 1.000 | 0.645 | 1 |
| 8 | SU11274 | 0.646 | 0.780 | 0.809 | 0.993 | 1.000 | 0.318 | 1 |
| 9 | H-89 | 0.616 | 0.633 | 0.667 | 0.872 | 1.000 | 0.893 | 1 |
| 10 | STO609 | 0.594 | 0.550 | 0.612 | 0.986 | 1.000 | 0.892 | 2 |
