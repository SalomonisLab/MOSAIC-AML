# MOSAIC-Rx drug prioritisation

*Predicted **ex-vivo** sensitivity from the BeatAML2 functional screen. This is a prioritisation for trial matching or laboratory validation, **not** a treatment recommendation.*

- inhibitors modelled: **118**, reported 118, abstained 0
- genes shared with the model: 14237
- distance from the BeatAML training distribution: 57th percentile
- differentiation axis (primitive - mature): 30th percentile of BeatAML
- cell states scored: 26

## approved in AML  (6 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Venetoclax | 0.782 | 0.962 | 0.966 | 0.667 | 1.000 | 0.434 | 2 |
| 2 | Azacytidine | 0.510 | 0.508 | 0.646 | 0.608 | 0.810 | 0.456 | 2 |
| 3 | Gilteritinib | 0.322 | 0.143 | 0.850 | 0.442 | 0.119 | 0.517 | 1 |
| 4 | Quizartinib (AC220) | 0.310 | 0.119 | 0.855 | 0.429 | 0.000 | 0.587 | 1 |
| 5 | Midostaurin | 0.283 | 0.181 | 0.803 | 0.236 | 0.000 | 0.587 | 2 |
| 6 | Cytarabine | 0.217 | 0.181 | 0.824 | 0.000 | 0.000 | 0.392 | 2 |

## approved, other indication (off-label in AML)  (36 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Afatinib (BIBW-2992) | 0.705 | 0.717 | 0.718 | 1.000 | 1.000 | 0.654 | 0 |
| 2 | Entrectinib | 0.662 | 0.691 | 0.771 | 0.931 | 1.000 | 0.535 | 1 |
| 3 | Neratinib (HKI-272) | 0.658 | 0.818 | 0.716 | 0.777 | 0.810 | 0.654 | 1 |
| 4 | Idelalisib | 0.631 | 0.793 | 0.759 | 0.846 | 0.714 | 0.599 | 0 |
| 5 | Lenalidomide | 0.628 | 0.697 | 0.726 | 0.667 | 1.000 | 0.588 | 1 |
| 6 | Pazopanib (GW786034) | 0.544 | 0.568 | 0.885 | 0.816 | 0.714 | 0.588 | 0 |
| 7 | Imatinib | 0.527 | 0.591 | 0.707 | 0.866 | 0.714 | 0.497 | 0 |
| 8 | Crizotinib (PF-2341066) | 0.474 | 0.784 | 0.802 | 0.730 | 0.278 | 0.429 | 1 |
| 9 | Selinexor | 0.436 | 0.526 | 0.550 | 0.489 | 0.810 | 0.455 | 2 |
| 10 | Lenvatinib | 0.371 | 0.457 | 0.861 | 0.638 | 0.000 | 0.555 | 0 |

## clinical-trial agent  (36 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Pelitinib (EKB-569) | 0.605 | 0.730 | 0.893 | 0.667 | 1.000 | 0.496 | 1 |
| 2 | Dovitinib (CHIR-258) | 0.584 | 0.818 | 0.917 | 0.747 | 0.714 | 0.505 | 1 |
| 3 | Elesclomol | 0.527 | 0.733 | 0.800 | 0.690 | 0.690 | 0.567 | 1 |
| 4 | Volasertib (BI-6727) | 0.466 | 0.506 | 0.668 | 0.886 | 1.000 | 0.158 | 0 |
| 5 | MGCD-265 | 0.461 | 0.861 | 0.843 | 0.546 | 0.373 | 0.432 | 2 |
| 6 | GDC-0941 | 0.382 | 0.449 | 0.825 | 0.707 | 0.413 | 0.411 | 1 |
| 7 | Canertinib (CI-1033) | 0.355 | 0.671 | 0.773 | 0.347 | 0.278 | 0.654 | 1 |
| 8 | GSK690693 | 0.332 | 0.551 | 0.658 | 0.774 | 0.278 | 0.476 | 0 |
| 9 | KW-2449 | 0.311 | 0.557 | 0.913 | 0.442 | 0.278 | 0.497 | 0 |
| 10 | MK-2206 | 0.292 | 0.477 | 0.767 | 0.529 | 0.000 | 0.476 | 0 |

## research-only compound  (40 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | AGI-6780 | 0.625 | 0.771 | 0.773 | 0.836 | 1.000 | 0.586 | 1 |
| 2 | Bay 11-7085 | 0.593 | 0.675 | 0.705 | 1.000 | 1.000 | 0.680 | 1 |
| 3 | AGI-5198 | 0.526 | 0.618 | 0.649 | 0.940 | 0.810 | 0.617 | 1 |
| 4 | JNJ-38877605 | 0.489 | 0.666 | 0.709 | 0.903 | 0.690 | 0.429 | 1 |
| 5 | ABT-737 | 0.475 | 0.759 | 0.763 | 0.571 | 0.690 | 0.451 | 1 |
| 6 | SR9011 | 0.458 | 0.655 | 0.721 | 0.747 | 0.714 | 0.608 | 1 |
| 7 | PI-103 | 0.446 | 0.705 | 0.757 | 0.764 | 0.595 | 0.342 | 1 |
| 8 | GSK-1838705A | 0.420 | 0.496 | 0.744 | 0.667 | 1.000 | 0.411 | 1 |
| 9 | DBZ | 0.408 | 0.565 | 0.667 | 0.553 | 0.881 | 0.477 | 1 |
| 10 | PHA-665752 | 0.405 | 0.810 | 0.798 | 0.471 | 0.373 | 0.429 | 1 |

## Combination hypotheses (complementary cell-state coverage)

*BeatAML2 measured single agents only; no synergy is claimed.*

| pair | pathways | coverage A | coverage B | union | gain |
|---|---|---|---|---|---|
| Axitinib (AG-013736) + Palbociclib | RTK / cell_cycle | 0.43 | 0.51 | 0.89 | +0.38 |
| CYT387 + Palbociclib | JAK_STAT / cell_cycle | 0.43 | 0.51 | 0.89 | +0.38 |
| Dasatinib + Palbociclib | RTK / cell_cycle | 0.43 | 0.51 | 0.89 | +0.38 |
| Erlotinib + Panobinostat | RTK / epigenetic | 0.41 | 0.38 | 0.78 | +0.38 |
| Ibrutinib (PCI-32765) + Palbociclib | immune_signalling / cell_cycle | 0.51 | 0.51 | 0.89 | +0.38 |
| Nilotinib + Palbociclib | RTK / cell_cycle | 0.43 | 0.51 | 0.89 | +0.38 |
| Palbociclib + Panobinostat | cell_cycle / epigenetic | 0.51 | 0.38 | 0.89 | +0.38 |
| CYT387 + Erlotinib | JAK_STAT / RTK | 0.43 | 0.41 | 0.78 | +0.35 |
