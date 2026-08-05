# MOSAIC-Rx drug prioritisation

*Predicted **ex-vivo** sensitivity from the BeatAML2 functional screen. This is a prioritisation for trial matching or laboratory validation, **not** a treatment recommendation.*

- inhibitors modelled: **118**, reported 118, abstained 0
- genes shared with the model: 14237
- distance from the BeatAML training distribution: 30th percentile
- differentiation axis (primitive - mature): 43th percentile of BeatAML
- cell states scored: 52

## approved in AML  (6 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Cytarabine | 0.835 | 0.837 | 0.824 | 0.939 | 0.991 | 0.715 | 1 |
| 2 | Quizartinib (AC220) | 0.735 | 0.632 | 0.855 | 0.875 | 0.985 | 0.651 | 0 |
| 3 | Azacytidine | 0.721 | 0.571 | 0.646 | 0.960 | 1.000 | 0.768 | 1 |
| 4 | Venetoclax | 0.715 | 0.675 | 0.966 | 0.867 | 1.000 | 0.550 | 1 |
| 5 | Gilteritinib | 0.493 | 0.557 | 0.850 | 0.514 | 0.477 | 0.548 | 1 |
| 6 | Midostaurin | 0.426 | 0.492 | 0.803 | 0.271 | 0.201 | 0.564 | 1 |

## approved, other indication (off-label in AML)  (36 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Ponatinib (AP24534) | 0.855 | 0.992 | 0.952 | 0.992 | 1.000 | 0.576 | 1 |
| 2 | Erlotinib | 0.839 | 0.881 | 0.834 | 0.958 | 1.000 | 0.864 | 0 |
| 3 | Sorafenib | 0.822 | 0.945 | 0.919 | 0.958 | 0.991 | 0.552 | 0 |
| 4 | Sunitinib | 0.814 | 0.858 | 0.891 | 0.937 | 0.991 | 0.759 | 1 |
| 5 | Cabozantinib | 0.811 | 0.888 | 0.901 | 0.896 | 1.000 | 0.685 | 0 |
| 6 | Crizotinib (PF-2341066) | 0.807 | 0.897 | 0.802 | 0.990 | 1.000 | 0.625 | 1 |
| 7 | Vargetef | 0.805 | 0.872 | 0.827 | 0.958 | 1.000 | 0.690 | 0 |
| 8 | Vandetanib (ZD6474) | 0.800 | 0.855 | 0.825 | 0.972 | 0.991 | 0.707 | 0 |
| 9 | Tivozanib (AV-951) | 0.796 | 0.918 | 0.941 | 0.903 | 0.888 | 0.667 | 0 |
| 10 | Afatinib (BIBW-2992) | 0.784 | 0.807 | 0.718 | 0.995 | 1.000 | 0.755 | 0 |

## clinical-trial agent  (36 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Pelitinib (EKB-569) | 0.835 | 0.958 | 0.893 | 0.964 | 1.000 | 0.864 | 0 |
| 2 | KW-2449 | 0.804 | 0.949 | 0.913 | 0.972 | 1.000 | 0.667 | 0 |
| 3 | Dovitinib (CHIR-258) | 0.771 | 0.863 | 0.917 | 0.928 | 0.994 | 0.731 | 1 |
| 4 | AZD1152-HQPA (AZD2811) | 0.766 | 0.842 | 0.824 | 0.919 | 0.917 | 0.968 | 1 |
| 5 | Alisertib (MLN8237) | 0.765 | 0.789 | 0.728 | 0.995 | 1.000 | 0.930 | 0 |
| 6 | Foretinib (XL880) | 0.761 | 0.896 | 0.908 | 0.904 | 0.991 | 0.619 | 0 |
| 7 | Canertinib (CI-1033) | 0.757 | 0.845 | 0.773 | 0.961 | 1.000 | 0.755 | 1 |
| 8 | OTX-015 | 0.741 | 0.849 | 0.779 | 0.995 | 1.000 | 0.624 | 0 |
| 9 | RAF265 (CHIR-265) | 0.737 | 0.866 | 0.831 | 0.938 | 1.000 | 0.586 | 1 |
| 10 | Volasertib (BI-6727) | 0.720 | 0.711 | 0.668 | 0.958 | 0.988 | 0.959 | 0 |

## research-only compound  (40 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | NVP-TAE684 | 0.823 | 0.964 | 0.907 | 0.966 | 0.991 | 0.864 | 0 |
| 2 | GW-2580 | 0.789 | 0.923 | 0.863 | 0.967 | 0.991 | 0.796 | 1 |
| 3 | PHA-665752 | 0.733 | 0.883 | 0.798 | 0.945 | 1.000 | 0.625 | 0 |
| 4 | JQ1 | 0.720 | 0.845 | 0.763 | 0.995 | 0.994 | 0.624 | 0 |
| 5 | PD173955 | 0.714 | 0.828 | 0.794 | 0.968 | 0.991 | 0.623 | 0 |
| 6 | AGI-6780 | 0.710 | 0.762 | 0.773 | 0.947 | 1.000 | 0.836 | 1 |
| 7 | A-674563 | 0.690 | 0.842 | 0.777 | 0.997 | 1.000 | 0.411 | 1 |
| 8 | KI20227 | 0.689 | 0.792 | 0.808 | 0.978 | 1.000 | 0.733 | 2 |
| 9 | PI-103 | 0.671 | 0.799 | 0.757 | 0.992 | 0.991 | 0.438 | 1 |
| 10 | JNJ-38877605 | 0.669 | 0.751 | 0.709 | 0.948 | 1.000 | 0.625 | 1 |
