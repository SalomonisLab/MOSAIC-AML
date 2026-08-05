# MOSAIC-Rx drug prioritisation

*Predicted **ex-vivo** sensitivity from the BeatAML2 functional screen. This is a prioritisation for trial matching or laboratory validation, **not** a treatment recommendation.*

- inhibitors modelled: **118**, reported 118, abstained 0
- genes shared with the model: 14237
- distance from the BeatAML training distribution: 83th percentile
- differentiation axis (primitive - mature): 33th percentile of BeatAML
- cell states scored: 41

## approved in AML  (6 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Venetoclax | 0.673 | 0.937 | 0.966 | 0.770 | 1.000 | 0.533 | 2 |
| 2 | Midostaurin | 0.259 | 0.355 | 0.803 | 0.363 | 0.231 | 0.446 | 0 |
| 3 | Cytarabine | 0.166 | 0.341 | 0.824 | 0.094 | 0.037 | 0.502 | 2 |
| 4 | Gilteritinib | 0.152 | 0.096 | 0.850 | 0.215 | 0.135 | 0.405 | 0 |
| 5 | Quizartinib (AC220) | 0.128 | 0.166 | 0.855 | 0.050 | 0.000 | 0.437 | 0 |
| 6 | Azacytidine | 0.127 | 0.395 | 0.646 | 0.016 | 0.000 | 0.451 | 2 |

## approved, other indication (off-label in AML)  (36 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Bortezomib (Velcade) | 0.574 | 0.685 | 0.667 | 0.974 | 1.000 | 0.657 | 1 |
| 2 | Neratinib (HKI-272) | 0.532 | 0.769 | 0.716 | 0.948 | 0.964 | 0.239 | 2 |
| 3 | Bosutinib (SKI-606) | 0.462 | 0.771 | 0.862 | 0.833 | 0.640 | 0.356 | 1 |
| 4 | Ponatinib (AP24534) | 0.440 | 0.548 | 0.952 | 0.648 | 0.920 | 0.464 | 0 |
| 5 | Tofacitinib (CP-690550) | 0.399 | 0.562 | 0.660 | 0.841 | 0.759 | 0.432 | 0 |
| 6 | Axitinib (AG-013736) | 0.394 | 0.684 | 0.862 | 0.798 | 0.522 | 0.376 | 0 |
| 7 | Entrectinib | 0.354 | 0.481 | 0.771 | 0.682 | 0.812 | 0.275 | 1 |
| 8 | Palbociclib | 0.350 | 0.629 | 0.782 | 0.616 | 0.657 | 0.547 | 0 |
| 9 | Lenalidomide | 0.334 | 0.530 | 0.726 | 0.520 | 0.657 | 0.559 | 1 |
| 10 | Trametinib (GSK1120212) | 0.317 | 0.409 | 0.903 | 0.636 | 0.768 | 0.165 | 0 |

## clinical-trial agent  (36 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | OTX-015 | 0.567 | 0.785 | 0.779 | 0.914 | 0.991 | 0.654 | 0 |
| 2 | Volasertib (BI-6727) | 0.555 | 0.727 | 0.668 | 0.992 | 0.991 | 0.703 | 0 |
| 3 | Cediranib (AZD2171) | 0.487 | 0.687 | 0.748 | 0.950 | 1.000 | 0.376 | 0 |
| 4 | Motesanib (AMG-706) | 0.470 | 0.648 | 0.824 | 0.973 | 0.965 | 0.345 | 0 |
| 5 | MK-2206 | 0.454 | 0.593 | 0.767 | 0.982 | 0.965 | 0.423 | 0 |
| 6 | GSK690693 | 0.425 | 0.556 | 0.658 | 0.974 | 0.955 | 0.423 | 0 |
| 7 | AZD1480 | 0.351 | 0.547 | 0.680 | 0.884 | 0.812 | 0.291 | 0 |
| 8 | Elesclomol | 0.309 | 0.780 | 0.800 | 0.728 | 0.692 | 0.017 | 2 |
| 9 | Masitinib (AB-1010) | 0.290 | 0.679 | 0.684 | 0.761 | 0.416 | 0.335 | 0 |
| 10 | Vatalanib (PTK787) | 0.275 | 0.618 | 0.727 | 0.771 | 0.416 | 0.376 | 0 |

## research-only compound  (40 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Nutlin 3a | 0.488 | 0.733 | 0.786 | 0.751 | 1.000 | 0.546 | 1 |
| 2 | JQ1 | 0.417 | 0.767 | 0.763 | 0.727 | 0.621 | 0.654 | 1 |
| 3 | JAK Inhibitor I | 0.357 | 0.460 | 0.700 | 0.914 | 0.823 | 0.450 | 0 |
| 4 | ABT-737 | 0.350 | 0.556 | 0.763 | 0.575 | 0.812 | 0.633 | 1 |
| 5 | BMS-345541 | 0.344 | 0.555 | 0.556 | 0.971 | 0.947 | 0.310 | 2 |
| 6 | PI-103 | 0.175 | 0.418 | 0.757 | 0.634 | 0.211 | 0.288 | 0 |
| 7 | DBZ | 0.172 | 0.501 | 0.667 | 0.419 | 0.605 | 0.420 | 2 |
| 8 | MLN8054 | 0.165 | 0.420 | 0.593 | 0.216 | 0.188 | 0.757 | 1 |
| 9 | PHT-427 | 0.161 | 0.516 | 0.598 | 0.746 | 0.371 | 0.466 | 2 |
| 10 | PRT062607 | 0.152 | 0.439 | 0.804 | 0.578 | 0.036 | 0.330 | 0 |

## Combination hypotheses (complementary cell-state coverage)

*BeatAML2 measured single agents only; no synergy is claimed.*

| pair | pathways | coverage A | coverage B | union | gain |
|---|---|---|---|---|---|
| CYT387 + Entrectinib | JAK_STAT / RTK | 0.56 | 0.56 | 0.88 | +0.32 |
| Afatinib (BIBW-2992) + Panobinostat | RTK / epigenetic | 0.26 | 0.24 | 0.49 | +0.24 |
| Afatinib (BIBW-2992) + Ibrutinib (PCI-32765) | RTK / immune_signalling | 0.26 | 0.22 | 0.48 | +0.22 |
| Entrectinib + Midostaurin | RTK / FLT3 | 0.56 | 0.62 | 0.83 | +0.22 |
| CYT387 + Palbociclib | JAK_STAT / cell_cycle | 0.56 | 0.81 | 1.00 | +0.20 |
| Cabozantinib + Gilteritinib | RTK / FLT3 | 0.32 | 0.34 | 0.54 | +0.20 |
| Crizotinib (PF-2341066) + Palbociclib | RTK / cell_cycle | 0.20 | 0.81 | 1.00 | +0.20 |
| Dasatinib + Palbociclib | RTK / cell_cycle | 0.22 | 0.81 | 1.00 | +0.20 |
