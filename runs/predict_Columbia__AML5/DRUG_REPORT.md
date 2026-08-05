# MOSAIC-Rx drug prioritisation

*Predicted **ex-vivo** sensitivity from the BeatAML2 functional screen. This is a prioritisation for trial matching or laboratory validation, **not** a treatment recommendation.*

- inhibitors modelled: **118**, reported 118, abstained 0
- genes shared with the model: 14237
- distance from the BeatAML training distribution: 37th percentile
- differentiation axis (primitive - mature): 64th percentile of BeatAML
- cell states scored: 22

## approved in AML  (6 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Venetoclax | 0.844 | 0.928 | 0.966 | 0.987 | 1.000 | 0.528 | 0 |
| 2 | Quizartinib (AC220) | 0.837 | 0.881 | 0.855 | 0.997 | 1.000 | 0.493 | 0 |
| 3 | Midostaurin | 0.723 | 0.662 | 0.803 | 0.922 | 0.997 | 0.481 | 1 |
| 4 | Gilteritinib | 0.688 | 0.524 | 0.850 | 0.920 | 0.997 | 0.586 | 1 |
| 5 | Azacytidine | 0.261 | 0.349 | 0.646 | 0.000 | 0.000 | 0.352 | 2 |
| 6 | Cytarabine | 0.244 | 0.302 | 0.824 | 0.005 | 0.000 | 0.217 | 2 |

## approved, other indication (off-label in AML)  (36 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Ponatinib (AP24534) | 0.835 | 0.981 | 0.952 | 0.980 | 1.000 | 0.501 | 2 |
| 2 | Cabozantinib | 0.822 | 0.944 | 0.901 | 0.933 | 1.000 | 0.597 | 0 |
| 3 | Sorafenib | 0.808 | 0.947 | 0.919 | 0.925 | 1.000 | 0.501 | 1 |
| 4 | Sunitinib | 0.805 | 0.923 | 0.891 | 0.933 | 1.000 | 0.545 | 0 |
| 5 | Erlotinib | 0.774 | 0.805 | 0.834 | 0.930 | 0.997 | 0.688 | 1 |
| 6 | Vargetef | 0.768 | 0.820 | 0.827 | 0.934 | 0.997 | 0.620 | 0 |
| 7 | Regorafenib (BAY 73-4506) | 0.734 | 0.777 | 0.834 | 0.933 | 1.000 | 0.519 | 0 |
| 8 | Gefitinib | 0.697 | 0.632 | 0.794 | 0.926 | 0.997 | 0.688 | 0 |
| 9 | Lenvatinib | 0.676 | 0.593 | 0.861 | 0.973 | 0.997 | 0.579 | 0 |
| 10 | Entrectinib | 0.670 | 0.687 | 0.771 | 0.747 | 1.000 | 0.614 | 1 |

## clinical-trial agent  (36 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Foretinib (XL880) | 0.791 | 0.939 | 0.908 | 0.939 | 1.000 | 0.648 | 1 |
| 2 | Dovitinib (CHIR-258) | 0.763 | 0.901 | 0.917 | 0.922 | 1.000 | 0.602 | 1 |
| 3 | KW-2449 | 0.756 | 0.927 | 0.913 | 0.925 | 1.000 | 0.490 | 1 |
| 4 | Crenolanib | 0.685 | 0.744 | 0.879 | 0.922 | 1.000 | 0.564 | 1 |
| 5 | Alisertib (MLN8237) | 0.617 | 0.724 | 0.728 | 0.995 | 0.997 | 0.203 | 1 |
| 6 | YM-155 | 0.614 | 0.608 | 0.754 | 0.921 | 1.000 | 0.530 | 1 |
| 7 | LY-333531 | 0.612 | 0.707 | 0.871 | 0.897 | 0.997 | 0.230 | 2 |
| 8 | Tozasertib (VX-680) | 0.606 | 0.665 | 0.728 | 0.936 | 0.997 | 0.347 | 0 |
| 9 | INK-128 | 0.583 | 0.578 | 0.793 | 0.988 | 1.000 | 0.319 | 1 |
| 10 | Linifanib (ABT-869) | 0.529 | 0.604 | 0.903 | 0.720 | 0.716 | 0.570 | 0 |

## research-only compound  (40 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | NVP-TAE684 | 0.719 | 0.801 | 0.907 | 0.918 | 0.997 | 0.688 | 0 |
| 2 | A-674563 | 0.629 | 0.732 | 0.777 | 0.980 | 1.000 | 0.347 | 0 |
| 3 | NF-kB Activation Inhibitor | 0.554 | 0.656 | 0.865 | 0.745 | 1.000 | 0.452 | 2 |
| 4 | GSK-1838705A | 0.544 | 0.511 | 0.744 | 0.892 | 0.997 | 0.699 | 2 |
| 5 | Bay 11-7085 | 0.538 | 0.536 | 0.705 | 0.983 | 1.000 | 0.545 | 1 |
| 6 | AGI-6780 | 0.530 | 0.475 | 0.773 | 0.752 | 0.997 | 0.616 | 1 |
| 7 | AGI-5198 | 0.525 | 0.474 | 0.649 | 0.829 | 0.997 | 0.590 | 1 |
| 8 | JAK Inhibitor I | 0.505 | 0.528 | 0.700 | 0.750 | 0.997 | 0.383 | 0 |
| 9 | KI20227 | 0.431 | 0.550 | 0.808 | 0.718 | 0.713 | 0.503 | 2 |
| 10 | KU-55933 | 0.429 | 0.640 | 0.737 | 0.777 | 0.713 | 0.236 | 3 |

## Combination hypotheses (complementary cell-state coverage)

*BeatAML2 measured single agents only; no synergy is claimed.*

| pair | pathways | coverage A | coverage B | union | gain |
|---|---|---|---|---|---|
| Ibrutinib (PCI-32765) + Imatinib | immune_signalling / RTK | 0.76 | 0.73 | 0.95 | +0.18 |
| Ibrutinib (PCI-32765) + Lenvatinib | immune_signalling / RTK | 0.76 | 0.40 | 0.95 | +0.18 |
| Afatinib (BIBW-2992) + Ibrutinib (PCI-32765) | RTK / immune_signalling | 0.16 | 0.76 | 0.92 | +0.16 |
| Dasatinib + Ibrutinib (PCI-32765) | RTK / immune_signalling | 0.16 | 0.76 | 0.92 | +0.16 |
