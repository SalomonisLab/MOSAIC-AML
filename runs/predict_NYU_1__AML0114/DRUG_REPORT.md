# MOSAIC-Rx drug prioritisation

*Predicted **ex-vivo** sensitivity from the BeatAML2 functional screen. This is a prioritisation for trial matching or laboratory validation, **not** a treatment recommendation.*

- inhibitors modelled: **118**, reported 118, abstained 0
- genes shared with the model: 14237
- distance from the BeatAML training distribution: 16th percentile
- differentiation axis (primitive - mature): 84th percentile of BeatAML
- cell states scored: 42

## approved in AML  (6 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Quizartinib (AC220) | 0.686 | 0.626 | 0.855 | 0.849 | 0.925 | 0.498 | 0 |
| 2 | Venetoclax | 0.609 | 0.419 | 0.966 | 0.744 | 0.929 | 0.657 | 1 |
| 3 | Midostaurin | 0.456 | 0.321 | 0.803 | 0.463 | 0.583 | 0.419 | 0 |
| 4 | Cytarabine | 0.420 | 0.347 | 0.824 | 0.239 | 0.254 | 0.711 | 1 |
| 5 | Gilteritinib | 0.392 | 0.365 | 0.850 | 0.218 | 0.286 | 0.483 | 1 |
| 6 | Azacytidine | 0.277 | 0.224 | 0.646 | 0.000 | 0.000 | 0.642 | 3 |

## approved, other indication (off-label in AML)  (36 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Ruxolitinib (INCB018424) | 0.723 | 0.765 | 0.713 | 0.939 | 1.000 | 0.575 | 0 |
| 2 | Gefitinib | 0.669 | 0.692 | 0.794 | 0.881 | 0.929 | 0.545 | 0 |
| 3 | CYT387 | 0.668 | 0.679 | 0.703 | 0.950 | 0.925 | 0.559 | 0 |
| 4 | Palbociclib | 0.652 | 0.652 | 0.782 | 0.802 | 0.925 | 0.646 | 0 |
| 5 | Tofacitinib (CP-690550) | 0.644 | 0.631 | 0.660 | 0.935 | 0.925 | 0.588 | 0 |
| 6 | Vandetanib (ZD6474) | 0.640 | 0.660 | 0.825 | 0.646 | 0.929 | 0.665 | 1 |
| 7 | Pazopanib (GW786034) | 0.636 | 0.682 | 0.885 | 0.620 | 0.929 | 0.563 | 1 |
| 8 | Erlotinib | 0.632 | 0.673 | 0.834 | 0.606 | 0.961 | 0.545 | 1 |
| 9 | Panobinostat | 0.628 | 0.634 | 0.860 | 0.838 | 0.859 | 0.553 | 0 |
| 10 | Sorafenib | 0.612 | 0.626 | 0.919 | 0.650 | 0.925 | 0.511 | 1 |

## clinical-trial agent  (36 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | RAF265 (CHIR-265) | 0.695 | 0.788 | 0.831 | 0.807 | 0.929 | 0.746 | 0 |
| 2 | YM-155 | 0.654 | 0.657 | 0.754 | 0.946 | 0.925 | 0.749 | 0 |
| 3 | CI-1040 (PD184352) | 0.651 | 0.687 | 0.819 | 0.949 | 0.925 | 0.606 | 0 |
| 4 | AT7519 | 0.632 | 0.715 | 0.787 | 0.867 | 0.925 | 0.541 | 1 |
| 5 | Saracatinib (AZD0530) | 0.602 | 0.691 | 0.872 | 0.680 | 0.929 | 0.555 | 1 |
| 6 | OTX-015 | 0.578 | 0.575 | 0.779 | 0.895 | 0.931 | 0.541 | 0 |
| 7 | Tozasertib (VX-680) | 0.565 | 0.660 | 0.728 | 0.592 | 0.929 | 0.577 | 0 |
| 8 | Volasertib (BI-6727) | 0.544 | 0.507 | 0.668 | 0.889 | 0.929 | 0.582 | 0 |
| 9 | AZD1480 | 0.532 | 0.497 | 0.680 | 0.852 | 0.925 | 0.534 | 0 |
| 10 | Pelitinib (EKB-569) | 0.525 | 0.548 | 0.893 | 0.544 | 0.929 | 0.545 | 0 |

## research-only compound  (40 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | NVP-TAE684 | 0.717 | 0.912 | 0.907 | 0.880 | 0.925 | 0.545 | 1 |
| 2 | ABT-737 | 0.629 | 0.680 | 0.763 | 0.804 | 0.964 | 0.751 | 1 |
| 3 | JNJ-7706621 | 0.605 | 0.659 | 0.663 | 0.952 | 0.925 | 0.598 | 0 |
| 4 | PD173955 | 0.568 | 0.673 | 0.794 | 0.642 | 0.929 | 0.562 | 1 |
| 5 | PRT062607 | 0.568 | 0.634 | 0.804 | 0.881 | 0.925 | 0.435 | 0 |
| 6 | GSK-1838705A | 0.560 | 0.741 | 0.744 | 0.752 | 0.929 | 0.474 | 2 |
| 7 | SU11274 | 0.554 | 0.768 | 0.809 | 0.801 | 0.929 | 0.625 | 0 |
| 8 | Doramapimod (BIRB 796) | 0.546 | 0.643 | 0.841 | 0.895 | 0.925 | 0.456 | 1 |
| 9 | Nutlin 3a | 0.533 | 0.449 | 0.786 | 0.766 | 0.925 | 0.744 | 1 |
| 10 | S31-201 | 0.510 | 0.586 | 0.673 | 0.721 | 0.929 | 0.637 | 1 |

## Combination hypotheses (complementary cell-state coverage)

*BeatAML2 measured single agents only; no synergy is claimed.*

| pair | pathways | coverage A | coverage B | union | gain |
|---|---|---|---|---|---|
| Midostaurin + Panobinostat | FLT3 / epigenetic | 0.48 | 0.70 | 0.98 | +0.28 |
| Afatinib (BIBW-2992) + Panobinostat | RTK / epigenetic | 0.66 | 0.70 | 0.96 | +0.26 |
| Cabozantinib + Panobinostat | RTK / epigenetic | 0.46 | 0.70 | 0.94 | +0.24 |
| Gilteritinib + Panobinostat | FLT3 / epigenetic | 0.33 | 0.70 | 0.93 | +0.23 |
| Bosutinib (SKI-606) + Panobinostat | RTK / epigenetic | 0.44 | 0.70 | 0.92 | +0.21 |
| Entrectinib + Panobinostat | RTK / epigenetic | 0.66 | 0.70 | 0.92 | +0.21 |
| Imatinib + Panobinostat | RTK / epigenetic | 0.44 | 0.70 | 0.92 | +0.21 |
| Cytarabine + Entrectinib | chemotherapy / RTK | 0.75 | 0.66 | 0.96 | +0.21 |
