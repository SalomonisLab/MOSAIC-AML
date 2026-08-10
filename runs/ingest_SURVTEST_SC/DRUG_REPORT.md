# COMPASS-AML drug prioritisation

*Predicted **ex-vivo** sensitivity from the BeatAML2 functional screen. This is a prioritisation for trial matching or laboratory validation, **not** a treatment recommendation.*

- inhibitors modelled: **118**, reported 118, abstained 0
- genes shared with the model: 13432
- distance from the BeatAML training distribution: 41th percentile
- differentiation axis (primitive - mature): 24th percentile of BeatAML
- cell states scored: 16

## approved in AML  (6 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Quizartinib (AC220) | 0.757 | 0.843 | 0.855 | 0.993 | 0.974 | 0.691 | 0 |
| 2 | Midostaurin | 0.701 | 0.815 | 0.803 | 0.957 | 0.844 | 0.693 | 0 |
| 3 | Gilteritinib | 0.347 | 0.557 | 0.850 | 0.556 | 0.029 | 0.688 | 0 |
| 4 | Cytarabine | 0.189 | 0.230 | 0.824 | 0.218 | 0.000 | 0.185 | 1 |
| 5 | Azacytidine | 0.186 | 0.242 | 0.646 | 0.000 | 0.000 | 0.427 | 1 |
| 6 | Venetoclax | 0.152 | 0.007 | 0.966 | 0.066 | 0.129 | 0.230 | 0 |

## approved, other indication (off-label in AML)  (36 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Trametinib (GSK1120212) | 0.868 | 0.950 | 0.903 | 1.000 | 1.000 | 0.791 | 0 |
| 2 | Tivozanib (AV-951) | 0.859 | 0.984 | 0.941 | 0.957 | 0.844 | 0.904 | 0 |
| 3 | Axitinib (AG-013736) | 0.857 | 0.914 | 0.862 | 0.993 | 0.974 | 0.872 | 0 |
| 4 | Regorafenib (BAY 73-4506) | 0.847 | 0.904 | 0.834 | 0.971 | 0.974 | 0.877 | 0 |
| 5 | Dasatinib | 0.841 | 0.975 | 0.931 | 0.957 | 0.844 | 0.822 | 0 |
| 6 | Vargetef | 0.839 | 0.866 | 0.827 | 0.993 | 0.974 | 0.899 | 0 |
| 7 | Pazopanib (GW786034) | 0.833 | 0.944 | 0.885 | 0.957 | 0.844 | 0.884 | 0 |
| 8 | Bosutinib (SKI-606) | 0.826 | 0.935 | 0.862 | 0.957 | 0.844 | 0.875 | 0 |
| 9 | Panobinostat | 0.821 | 0.918 | 0.860 | 1.000 | 1.000 | 0.618 | 0 |
| 10 | Lenvatinib | 0.821 | 0.915 | 0.861 | 0.957 | 0.844 | 0.902 | 0 |

## clinical-trial agent  (36 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Saracatinib (AZD0530) | 0.833 | 0.943 | 0.872 | 0.993 | 0.974 | 0.902 | 0 |
| 2 | Entospletinib (GS-9973) | 0.801 | 0.961 | 0.872 | 0.993 | 0.974 | 0.679 | 0 |
| 3 | CI-1040 (PD184352) | 0.793 | 0.888 | 0.819 | 1.000 | 1.000 | 0.791 | 0 |
| 4 | LY-333531 | 0.752 | 0.861 | 0.871 | 0.993 | 0.974 | 0.628 | 0 |
| 5 | Motesanib (AMG-706) | 0.742 | 0.847 | 0.824 | 0.957 | 0.844 | 0.878 | 0 |
| 6 | Cediranib (AZD2171) | 0.733 | 0.825 | 0.748 | 0.964 | 0.871 | 0.872 | 0 |
| 7 | OTX-015 | 0.717 | 0.778 | 0.779 | 0.978 | 0.971 | 0.694 | 0 |
| 8 | Linifanib (ABT-869) | 0.670 | 0.922 | 0.903 | 0.957 | 0.844 | 0.732 | 0 |
| 9 | Crenolanib | 0.655 | 0.816 | 0.879 | 0.993 | 0.974 | 0.686 | 0 |
| 10 | AZD1480 | 0.632 | 0.684 | 0.680 | 0.933 | 0.816 | 0.778 | 0 |

## research-only compound  (40 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | GW-2580 | 0.764 | 0.923 | 0.863 | 0.957 | 0.844 | 0.884 | 0 |
| 2 | JNJ-28312141 | 0.750 | 0.899 | 0.842 | 0.993 | 0.974 | 0.839 | 1 |
| 3 | Doramapimod (BIRB 796) | 0.727 | 0.909 | 0.841 | 0.993 | 0.974 | 0.687 | 1 |
| 4 | PD173955 | 0.709 | 0.815 | 0.794 | 0.941 | 0.844 | 0.875 | 0 |
| 5 | PRT062607 | 0.700 | 0.871 | 0.804 | 0.957 | 0.844 | 0.679 | 0 |
| 6 | SU11274 | 0.689 | 0.678 | 0.809 | 0.971 | 0.974 | 0.869 | 0 |
| 7 | KI20227 | 0.675 | 0.831 | 0.808 | 0.941 | 0.844 | 0.835 | 1 |
| 8 | S31-201 | 0.585 | 0.635 | 0.673 | 0.957 | 0.844 | 0.858 | 1 |
| 9 | H-89 | 0.573 | 0.675 | 0.667 | 0.949 | 0.816 | 0.735 | 1 |
| 10 | BMS-345541 | 0.570 | 0.557 | 0.556 | 0.992 | 0.971 | 0.793 | 2 |

## Combination hypotheses (complementary cell-state coverage)

*BeatAML2 measured single agents only; no synergy is claimed.*

| pair | pathways | coverage A | coverage B | union | gain |
|---|---|---|---|---|---|
| Erlotinib + Palbociclib | RTK / cell_cycle | 0.56 | 0.49 | 1.00 | +0.44 |
| Cabozantinib + Palbociclib | RTK / cell_cycle | 0.70 | 0.49 | 1.00 | +0.30 |
| Idelalisib + Palbociclib | PI3K_AKT_mTOR / cell_cycle | 0.73 | 0.49 | 1.00 | +0.27 |
