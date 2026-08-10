# COMPASS-AML drug prioritisation

*Predicted **ex-vivo** sensitivity from the BeatAML2 functional screen. This is a prioritisation for trial matching or laboratory validation, **not** a treatment recommendation.*

- inhibitors modelled: **118**, reported 117, abstained 1
- genes shared with the model: 14237
- distance from the BeatAML training distribution: 33th percentile
- differentiation axis (primitive - mature): 28th percentile of BeatAML
- cell states scored: 32

## approved in AML  (6 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Quizartinib (AC220) | 0.862 | 0.860 | 0.855 | 0.994 | 1.000 | 0.653 | 0 |
| 2 | Midostaurin | 0.844 | 0.820 | 0.803 | 0.986 | 1.000 | 0.659 | 0 |
| 3 | Gilteritinib | 0.457 | 0.632 | 0.850 | 0.644 | 0.036 | 0.646 | 1 |
| 4 | Cytarabine | 0.191 | 0.236 | 0.824 | 0.159 | 0.000 | 0.243 | 1 |
| 5 | Azacytidine | 0.134 | 0.237 | 0.646 | 0.000 | 0.000 | 0.244 | 2 |
| 6 | Venetoclax | 0.132 | 0.013 | 0.966 | 0.027 | 0.000 | 0.211 | 1 |

## approved, other indication (off-label in AML)  (35 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Pazopanib (GW786034) | 0.859 | 0.944 | 0.885 | 0.987 | 0.959 | 0.826 | 0 |
| 2 | Axitinib (AG-013736) | 0.852 | 0.914 | 0.862 | 0.998 | 1.000 | 0.812 | 0 |
| 3 | Tivozanib (AV-951) | 0.851 | 0.984 | 0.941 | 0.954 | 0.836 | 0.879 | 0 |
| 4 | Trametinib (GSK1120212) | 0.847 | 0.946 | 0.903 | 0.987 | 0.994 | 0.658 | 0 |
| 5 | Lenvatinib | 0.838 | 0.913 | 0.861 | 0.987 | 0.959 | 0.849 | 0 |
| 6 | Vargetef | 0.828 | 0.878 | 0.827 | 0.998 | 1.000 | 0.841 | 0 |
| 7 | Dasatinib | 0.806 | 0.974 | 0.931 | 0.965 | 0.877 | 0.697 | 0 |
| 8 | Regorafenib (BAY 73-4506) | 0.800 | 0.892 | 0.834 | 0.954 | 0.836 | 0.814 | 0 |
| 9 | Gefitinib | 0.796 | 0.812 | 0.794 | 0.987 | 0.959 | 0.836 | 0 |
| 10 | Bosutinib (SKI-606) | 0.795 | 0.935 | 0.862 | 0.954 | 0.836 | 0.781 | 0 |

## clinical-trial agent  (36 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Saracatinib (AZD0530) | 0.814 | 0.942 | 0.872 | 0.998 | 1.000 | 0.823 | 0 |
| 2 | CI-1040 (PD184352) | 0.776 | 0.885 | 0.819 | 0.993 | 0.994 | 0.658 | 0 |
| 3 | Entospletinib (GS-9973) | 0.776 | 0.959 | 0.872 | 0.998 | 1.000 | 0.660 | 0 |
| 4 | Linifanib (ABT-869) | 0.770 | 0.922 | 0.903 | 0.965 | 0.877 | 0.753 | 0 |
| 5 | Crenolanib | 0.769 | 0.849 | 0.879 | 0.994 | 1.000 | 0.650 | 0 |
| 6 | Motesanib (AMG-706) | 0.739 | 0.841 | 0.824 | 0.954 | 0.836 | 0.817 | 0 |
| 7 | Cediranib (AZD2171) | 0.736 | 0.820 | 0.748 | 0.965 | 0.877 | 0.812 | 0 |
| 8 | LY-333531 | 0.721 | 0.866 | 0.871 | 0.974 | 0.959 | 0.627 | 0 |
| 9 | GDC-0941 | 0.711 | 0.847 | 0.825 | 0.945 | 0.800 | 0.706 | 1 |
| 10 | YM-155 | 0.680 | 0.812 | 0.754 | 0.990 | 0.964 | 0.371 | 1 |

## research-only compound  (40 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | GW-2580 | 0.739 | 0.920 | 0.863 | 0.954 | 0.836 | 0.780 | 1 |
| 2 | JNJ-28312141 | 0.733 | 0.895 | 0.842 | 0.998 | 1.000 | 0.718 | 1 |
| 3 | Doramapimod (BIRB 796) | 0.717 | 0.907 | 0.841 | 0.985 | 0.953 | 0.690 | 1 |
| 4 | NVP-TAE684 | 0.699 | 0.835 | 0.907 | 0.927 | 0.806 | 0.836 | 1 |
| 5 | PRT062607 | 0.669 | 0.870 | 0.804 | 0.954 | 0.836 | 0.660 | 0 |
| 6 | PD173955 | 0.663 | 0.807 | 0.794 | 0.945 | 0.800 | 0.781 | 0 |
| 7 | KI20227 | 0.654 | 0.823 | 0.808 | 0.941 | 0.806 | 0.759 | 1 |
| 8 | SU11274 | 0.648 | 0.693 | 0.809 | 0.943 | 0.830 | 0.768 | 0 |
| 9 | JNJ-7706621 | 0.603 | 0.700 | 0.663 | 0.998 | 1.000 | 0.287 | 1 |
| 10 | H-89 | 0.558 | 0.682 | 0.667 | 0.946 | 0.830 | 0.723 | 1 |

## Combination hypotheses (complementary cell-state coverage)

*BeatAML2 measured single agents only; no synergy is claimed.*

| pair | pathways | coverage A | coverage B | union | gain |
|---|---|---|---|---|---|
| Crizotinib (PF-2341066) + Palbociclib | RTK / cell_cycle | 0.73 | 0.52 | 1.00 | +0.27 |
| Idelalisib + Palbociclib | PI3K_AKT_mTOR / cell_cycle | 0.73 | 0.52 | 1.00 | +0.27 |
