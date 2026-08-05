# MOSAIC-Rx drug prioritisation

*Predicted **ex-vivo** sensitivity from the BeatAML2 functional screen. This is a prioritisation for trial matching or laboratory validation, **not** a treatment recommendation.*

- inhibitors modelled: **118**, reported 118, abstained 0
- genes shared with the model: 14237
- distance from the BeatAML training distribution: 67th percentile
- differentiation axis (primitive - mature): 83th percentile of BeatAML
- cell states scored: 23

## approved in AML  (6 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Venetoclax | 0.673 | 0.945 | 0.966 | 0.887 | 1.000 | 0.167 | 2 |
| 2 | Quizartinib (AC220) | 0.615 | 0.521 | 0.855 | 0.898 | 1.000 | 0.573 | 1 |
| 3 | Midostaurin | 0.547 | 0.434 | 0.803 | 0.706 | 0.921 | 0.632 | 1 |
| 4 | Cytarabine | 0.543 | 0.644 | 0.824 | 0.857 | 0.981 | 0.212 | 2 |
| 5 | Gilteritinib | 0.224 | 0.276 | 0.850 | 0.000 | 0.000 | 0.520 | 1 |
| 6 | Azacytidine | 0.052 | 0.230 | 0.646 | 0.000 | 0.000 | 0.042 | 2 |

## approved, other indication (off-label in AML)  (36 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Lenvatinib | 0.772 | 0.902 | 0.861 | 0.985 | 1.000 | 0.806 | 0 |
| 2 | Gefitinib | 0.765 | 0.836 | 0.794 | 0.985 | 1.000 | 0.941 | 0 |
| 3 | Regorafenib (BAY 73-4506) | 0.762 | 0.882 | 0.834 | 0.938 | 1.000 | 0.834 | 1 |
| 4 | Panobinostat | 0.721 | 0.920 | 0.860 | 1.000 | 1.000 | 0.417 | 1 |
| 5 | Selumetinib (AZD6244) | 0.669 | 0.705 | 0.886 | 0.977 | 1.000 | 0.643 | 0 |
| 6 | Ibrutinib (PCI-32765) | 0.667 | 0.829 | 0.825 | 0.930 | 1.000 | 0.415 | 1 |
| 7 | Ruxolitinib (INCB018424) | 0.660 | 0.798 | 0.713 | 0.993 | 1.000 | 0.476 | 0 |
| 8 | Vandetanib (ZD6474) | 0.659 | 0.742 | 0.825 | 0.905 | 0.940 | 0.700 | 1 |
| 9 | Entrectinib | 0.649 | 0.663 | 0.771 | 0.922 | 1.000 | 0.797 | 1 |
| 10 | Nilotinib | 0.627 | 0.648 | 0.690 | 0.890 | 0.940 | 0.843 | 0 |

## clinical-trial agent  (36 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Saracatinib (AZD0530) | 0.729 | 0.876 | 0.872 | 0.985 | 1.000 | 0.824 | 0 |
| 2 | RAF265 (CHIR-265) | 0.718 | 0.871 | 0.831 | 0.978 | 1.000 | 0.806 | 1 |
| 3 | Entospletinib (GS-9973) | 0.682 | 0.949 | 0.872 | 0.973 | 1.000 | 0.383 | 2 |
| 4 | Cediranib (AZD2171) | 0.620 | 0.666 | 0.748 | 0.953 | 1.000 | 0.799 | 0 |
| 5 | CI-1040 (PD184352) | 0.618 | 0.700 | 0.819 | 0.977 | 1.000 | 0.643 | 0 |
| 6 | LY-333531 | 0.617 | 0.843 | 0.871 | 0.942 | 1.000 | 0.277 | 2 |
| 7 | OTX-015 | 0.611 | 0.762 | 0.779 | 0.985 | 1.000 | 0.492 | 0 |
| 8 | Pelitinib (EKB-569) | 0.603 | 0.637 | 0.893 | 0.788 | 0.940 | 0.941 | 1 |
| 9 | GDC-0941 | 0.603 | 0.674 | 0.825 | 0.927 | 0.981 | 0.699 | 2 |
| 10 | Vatalanib (PTK787) | 0.600 | 0.713 | 0.727 | 0.872 | 0.921 | 0.799 | 1 |

## research-only compound  (40 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | SU11274 | 0.694 | 0.875 | 0.809 | 0.993 | 1.000 | 0.752 | 0 |
| 2 | Doramapimod (BIRB 796) | 0.595 | 0.792 | 0.841 | 0.984 | 1.000 | 0.566 | 1 |
| 3 | JNJ-38877605 | 0.594 | 0.704 | 0.709 | 0.911 | 0.981 | 0.752 | 1 |
| 4 | Nutlin 3a | 0.591 | 0.756 | 0.786 | 0.972 | 1.000 | 0.486 | 1 |
| 5 | PRT062607 | 0.563 | 0.714 | 0.804 | 0.978 | 1.000 | 0.383 | 0 |
| 6 | H-89 | 0.550 | 0.730 | 0.667 | 1.000 | 1.000 | 0.544 | 1 |
| 7 | Artemisinin | 0.546 | 0.593 | 0.665 | 0.987 | 1.000 | — | 1 |
| 8 | TG100-115 | 0.540 | 0.620 | 0.689 | 0.931 | 0.981 | 0.647 | 0 |
| 9 | S31-201 | 0.516 | 0.646 | 0.673 | 0.970 | 1.000 | 0.578 | 1 |
| 10 | JNJ-7706621 | 0.513 | 0.723 | 0.663 | 0.993 | 1.000 | 0.155 | 1 |
