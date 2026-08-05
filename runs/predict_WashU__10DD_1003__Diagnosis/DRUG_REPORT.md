# MOSAIC-Rx drug prioritisation

*Predicted **ex-vivo** sensitivity from the BeatAML2 functional screen. This is a prioritisation for trial matching or laboratory validation, **not** a treatment recommendation.*

- inhibitors modelled: **118**, reported 117, abstained 1
- genes shared with the model: 14237
- distance from the BeatAML training distribution: 52th percentile
- differentiation axis (primitive - mature): 26th percentile of BeatAML
- cell states scored: 27

## approved in AML  (6 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Cytarabine | 0.609 | 0.820 | 0.824 | 0.713 | 0.578 | 0.304 | 3 |
| 2 | Azacytidine | 0.544 | 0.511 | 0.646 | 0.585 | 1.000 | 0.291 | 2 |
| 3 | Venetoclax | 0.479 | 0.058 | 0.966 | 0.629 | 1.000 | 0.361 | 1 |
| 4 | Midostaurin | 0.351 | 0.311 | 0.803 | 0.384 | 0.023 | 0.519 | 0 |
| 5 | Quizartinib (AC220) | 0.268 | 0.151 | 0.855 | 0.042 | 0.063 | 0.518 | 0 |
| 6 | Gilteritinib | 0.247 | 0.121 | 0.850 | 0.087 | 0.023 | 0.437 | 0 |

## approved, other indication (off-label in AML)  (36 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Idelalisib | 0.692 | 0.809 | 0.759 | 0.984 | 0.973 | 0.293 | 1 |
| 2 | Entrectinib | 0.616 | 0.623 | 0.771 | 0.909 | 0.977 | 0.367 | 1 |
| 3 | Afatinib (BIBW-2992) | 0.546 | 0.534 | 0.718 | 0.619 | 0.983 | 0.482 | 0 |
| 4 | Ponatinib (AP24534) | 0.524 | 0.573 | 0.952 | 0.709 | 0.694 | 0.470 | 1 |
| 5 | Erlotinib | 0.489 | 0.383 | 0.834 | 0.652 | 1.000 | 0.266 | 0 |
| 6 | Crizotinib (PF-2341066) | 0.478 | 0.644 | 0.802 | 0.433 | 0.738 | 0.293 | 1 |
| 7 | Lenalidomide | 0.464 | 0.636 | 0.726 | 0.629 | 1.000 | 0.253 | 1 |
| 8 | Cabozantinib | 0.413 | 0.362 | 0.901 | 0.624 | 0.500 | 0.320 | 0 |
| 9 | Regorafenib (BAY 73-4506) | 0.400 | 0.503 | 0.834 | 0.621 | 0.494 | 0.353 | 0 |
| 10 | Trametinib (GSK1120212) | 0.381 | 0.867 | 0.903 | 0.447 | 0.056 | 0.177 | 1 |

## clinical-trial agent  (36 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | GSK690693 | 0.537 | 0.547 | 0.658 | 0.941 | 0.989 | 0.326 | 0 |
| 2 | MGCD-265 | 0.528 | 0.742 | 0.843 | 0.779 | 0.694 | 0.301 | 2 |
| 3 | Volasertib (BI-6727) | 0.525 | 0.580 | 0.668 | 0.923 | 0.922 | 0.314 | 0 |
| 4 | Canertinib (CI-1033) | 0.473 | 0.433 | 0.773 | 0.579 | 0.989 | 0.482 | 1 |
| 5 | Pelitinib (EKB-569) | 0.412 | 0.287 | 0.893 | 0.579 | 0.989 | 0.266 | 0 |
| 6 | KW-2449 | 0.369 | 0.302 | 0.913 | 0.408 | 0.694 | 0.413 | 1 |
| 7 | OTX-015 | 0.310 | 0.744 | 0.779 | 0.384 | 0.023 | 0.475 | 0 |
| 8 | Motesanib (AMG-706) | 0.302 | 0.729 | 0.824 | 0.409 | 0.056 | 0.387 | 0 |
| 9 | Elesclomol | 0.299 | 0.471 | 0.800 | 0.313 | 0.069 | 0.576 | 1 |
| 10 | Dovitinib (CHIR-258) | 0.291 | 0.206 | 0.917 | 0.300 | 0.508 | 0.374 | 1 |

## research-only compound  (39 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | GW-2580 | 0.653 | 0.909 | 0.863 | 0.808 | 0.744 | 0.612 | 1 |
| 2 | Nutlin 3a | 0.561 | 0.625 | 0.786 | 0.932 | 0.986 | 0.338 | 2 |
| 3 | NVP-ADW742 | 0.555 | 0.670 | 0.723 | 0.967 | 1.000 | 0.342 | 1 |
| 4 | A-674563 | 0.553 | 0.668 | 0.777 | 0.881 | 0.973 | 0.259 | 1 |
| 5 | JQ1 | 0.532 | 0.760 | 0.763 | 0.747 | 0.694 | 0.475 | 1 |
| 6 | GSK-1904529A | 0.521 | 0.624 | 0.696 | 0.984 | 0.973 | 0.300 | 2 |
| 7 | PHA-665752 | 0.485 | 0.594 | 0.798 | 0.585 | 1.000 | 0.293 | 1 |
| 8 | STO609 | 0.473 | 0.511 | 0.612 | 0.912 | 1.000 | 0.397 | 2 |
| 9 | DBZ | 0.470 | 0.627 | 0.667 | 0.657 | 1.000 | 0.317 | 1 |
| 10 | PI-103 | 0.423 | 0.612 | 0.757 | 0.752 | 0.686 | 0.199 | 1 |

## Combination hypotheses (complementary cell-state coverage)

*BeatAML2 measured single agents only; no synergy is claimed.*

| pair | pathways | coverage A | coverage B | union | gain |
|---|---|---|---|---|---|
| Afatinib (BIBW-2992) + Bortezomib (Velcade) | RTK / proteostasis | 0.63 | 0.37 | 1.00 | +0.37 |
| Afatinib (BIBW-2992) + CYT387 | RTK / JAK_STAT | 0.63 | 0.37 | 1.00 | +0.37 |
| Afatinib (BIBW-2992) + Ibrutinib (PCI-32765) | RTK / immune_signalling | 0.63 | 0.37 | 1.00 | +0.37 |
| Afatinib (BIBW-2992) + Midostaurin | RTK / FLT3 | 0.63 | 0.37 | 1.00 | +0.37 |
| Afatinib (BIBW-2992) + Panobinostat | RTK / epigenetic | 0.63 | 0.37 | 1.00 | +0.37 |
| Bortezomib (Velcade) + Erlotinib | proteostasis / RTK | 0.37 | 0.66 | 1.00 | +0.34 |
| CYT387 + Erlotinib | JAK_STAT / RTK | 0.37 | 0.66 | 1.00 | +0.34 |
| Erlotinib + Ibrutinib (PCI-32765) | RTK / immune_signalling | 0.66 | 0.37 | 1.00 | +0.34 |
