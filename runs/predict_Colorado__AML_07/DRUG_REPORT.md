# MOSAIC-Rx drug prioritisation

*Predicted **ex-vivo** sensitivity from the BeatAML2 functional screen. This is a prioritisation for trial matching or laboratory validation, **not** a treatment recommendation.*

- inhibitors modelled: **118**, reported 118, abstained 0
- genes shared with the model: 14237
- distance from the BeatAML training distribution: 19th percentile
- differentiation axis (primitive - mature): 56th percentile of BeatAML
- cell states scored: 32

## approved in AML  (6 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Quizartinib (AC220) | 0.605 | 0.474 | 0.855 | 0.646 | 0.902 | 0.533 | 0 |
| 2 | Cytarabine | 0.547 | 0.583 | 0.824 | 0.543 | 0.698 | 0.591 | 1 |
| 3 | Venetoclax | 0.535 | 0.107 | 0.966 | 0.693 | 1.000 | 0.676 | 2 |
| 4 | Azacytidine | 0.440 | 0.430 | 0.646 | 0.359 | 0.310 | 0.651 | 1 |
| 5 | Midostaurin | 0.429 | 0.492 | 0.803 | 0.440 | 0.101 | 0.497 | 0 |
| 6 | Gilteritinib | 0.317 | 0.264 | 0.850 | 0.057 | 0.003 | 0.624 | 1 |

## approved, other indication (off-label in AML)  (36 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Ibrutinib (PCI-32765) | 0.824 | 0.872 | 0.825 | 1.000 | 1.000 | 0.748 | 0 |
| 2 | Entrectinib | 0.733 | 0.782 | 0.771 | 1.000 | 1.000 | 0.486 | 1 |
| 3 | Ruxolitinib (INCB018424) | 0.720 | 0.750 | 0.713 | 0.998 | 0.997 | 0.508 | 0 |
| 4 | Lenvatinib | 0.717 | 0.822 | 0.861 | 0.906 | 0.902 | 0.462 | 0 |
| 5 | CYT387 | 0.713 | 0.727 | 0.703 | 0.986 | 1.000 | 0.539 | 0 |
| 6 | Crizotinib (PF-2341066) | 0.712 | 0.840 | 0.802 | 0.807 | 1.000 | 0.372 | 1 |
| 7 | Gefitinib | 0.711 | 0.752 | 0.794 | 0.908 | 1.000 | 0.488 | 0 |
| 8 | Erlotinib | 0.697 | 0.791 | 0.834 | 0.700 | 1.000 | 0.488 | 0 |
| 9 | Sunitinib | 0.689 | 0.737 | 0.891 | 0.700 | 1.000 | 0.536 | 0 |
| 10 | Imatinib | 0.683 | 0.683 | 0.707 | 0.915 | 1.000 | 0.528 | 0 |

## clinical-trial agent  (36 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | LY-333531 | 0.763 | 0.830 | 0.871 | 0.906 | 0.997 | 0.823 | 1 |
| 2 | Pelitinib (EKB-569) | 0.708 | 0.911 | 0.893 | 0.709 | 1.000 | 0.488 | 1 |
| 3 | Vatalanib (PTK787) | 0.665 | 0.733 | 0.727 | 1.000 | 1.000 | 0.457 | 0 |
| 4 | OTX-015 | 0.664 | 0.811 | 0.779 | 0.858 | 0.788 | 0.715 | 0 |
| 5 | AZD1480 | 0.655 | 0.660 | 0.680 | 0.988 | 0.997 | 0.626 | 1 |
| 6 | RAF265 (CHIR-265) | 0.655 | 0.800 | 0.831 | 0.706 | 1.000 | 0.473 | 1 |
| 7 | Masitinib (AB-1010) | 0.633 | 0.655 | 0.684 | 0.974 | 1.000 | 0.516 | 0 |
| 8 | Saracatinib (AZD0530) | 0.626 | 0.682 | 0.872 | 0.877 | 0.832 | 0.641 | 0 |
| 9 | Volasertib (BI-6727) | 0.582 | 0.573 | 0.668 | 0.809 | 1.000 | 0.603 | 0 |
| 10 | Entospletinib (GS-9973) | 0.567 | 0.617 | 0.872 | 0.777 | 0.654 | 0.836 | 1 |

## research-only compound  (40 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | PRT062607 | 0.757 | 0.827 | 0.804 | 0.995 | 0.997 | 0.836 | 0 |
| 2 | JQ1 | 0.722 | 0.799 | 0.763 | 0.998 | 1.000 | 0.715 | 0 |
| 3 | GW-2580 | 0.700 | 0.766 | 0.863 | 0.887 | 1.000 | 0.710 | 0 |
| 4 | NVP-TAE684 | 0.677 | 0.867 | 0.907 | 0.709 | 1.000 | 0.488 | 0 |
| 5 | JAK Inhibitor I | 0.661 | 0.736 | 0.700 | 0.997 | 1.000 | 0.542 | 0 |
| 6 | Doramapimod (BIRB 796) | 0.658 | 0.788 | 0.841 | 0.957 | 1.000 | 0.550 | 1 |
| 7 | NF-kB Activation Inhibitor | 0.647 | 0.906 | 0.865 | 0.802 | 1.000 | 0.342 | 3 |
| 8 | PD173955 | 0.640 | 0.646 | 0.794 | 0.936 | 1.000 | 0.639 | 0 |
| 9 | GSK-1838705A | 0.606 | 0.761 | 0.744 | 0.701 | 1.000 | 0.661 | 2 |
| 10 | AGI-5198 | 0.593 | 0.574 | 0.649 | 0.827 | 1.000 | 0.759 | 1 |

## Combination hypotheses (complementary cell-state coverage)

*BeatAML2 measured single agents only; no synergy is claimed.*

| pair | pathways | coverage A | coverage B | union | gain |
|---|---|---|---|---|---|
| Afatinib (BIBW-2992) + Bortezomib (Velcade) | RTK / proteostasis | 0.70 | 0.30 | 1.00 | +0.30 |
| Afatinib (BIBW-2992) + Panobinostat | RTK / epigenetic | 0.70 | 0.30 | 1.00 | +0.30 |
| Axitinib (AG-013736) + Cytarabine | RTK / chemotherapy | 0.35 | 0.70 | 1.00 | +0.30 |
| Axitinib (AG-013736) + Lenalidomide | RTK / epigenetic | 0.35 | 0.70 | 1.00 | +0.30 |
| Axitinib (AG-013736) + Palbociclib | RTK / cell_cycle | 0.35 | 0.70 | 1.00 | +0.30 |
| Bortezomib (Velcade) + Cytarabine | proteostasis / chemotherapy | 0.30 | 0.70 | 1.00 | +0.30 |
| Bortezomib (Velcade) + Lapatinib | proteostasis / RTK | 0.30 | 0.70 | 1.00 | +0.30 |
| Bortezomib (Velcade) + Lenalidomide | proteostasis / epigenetic | 0.30 | 0.70 | 1.00 | +0.30 |
