# MOSAIC-Rx drug prioritisation

*Predicted **ex-vivo** sensitivity from the BeatAML2 functional screen. This is a prioritisation for trial matching or laboratory validation, **not** a treatment recommendation.*

- inhibitors modelled: **118**, reported 118, abstained 0
- genes shared with the model: 14237
- distance from the BeatAML training distribution: 16th percentile
- differentiation axis (primitive - mature): 5th percentile of BeatAML
- cell states scored: 38

## approved in AML  (6 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Cytarabine | 0.419 | 0.583 | 0.824 | 0.206 | 0.746 | 0.228 | 1 |
| 2 | Quizartinib (AC220) | 0.401 | 0.313 | 0.855 | 0.020 | 0.643 | 0.437 | 0 |
| 3 | Venetoclax | 0.317 | 0.003 | 0.966 | 0.031 | 1.000 | 0.185 | 2 |
| 4 | Midostaurin | 0.293 | 0.291 | 0.803 | 0.048 | 0.032 | 0.454 | 1 |
| 5 | Azacytidine | 0.257 | 0.330 | 0.646 | 0.037 | 0.203 | 0.317 | 2 |
| 6 | Gilteritinib | 0.252 | 0.152 | 0.850 | 0.003 | 0.015 | 0.457 | 1 |

## approved, other indication (off-label in AML)  (36 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Ibrutinib (PCI-32765) | 0.772 | 0.839 | 0.825 | 0.984 | 0.934 | 0.655 | 0 |
| 2 | Lenvatinib | 0.767 | 0.894 | 0.861 | 0.984 | 0.956 | 0.434 | 1 |
| 3 | Entrectinib | 0.689 | 0.770 | 0.771 | 0.982 | 0.956 | 0.364 | 2 |
| 4 | Ruxolitinib (INCB018424) | 0.680 | 0.765 | 0.713 | 0.995 | 0.847 | 0.479 | 0 |
| 5 | CYT387 | 0.660 | 0.736 | 0.703 | 0.993 | 0.800 | 0.513 | 0 |
| 6 | Rapamycin | 0.627 | 0.952 | 0.922 | 0.976 | 0.254 | 0.518 | 2 |
| 7 | Dasatinib | 0.626 | 0.954 | 0.931 | 0.956 | 0.254 | 0.503 | 1 |
| 8 | Ponatinib (AP24534) | 0.601 | 0.706 | 0.952 | 0.693 | 0.746 | 0.458 | 0 |
| 9 | Bosutinib (SKI-606) | 0.582 | 0.616 | 0.862 | 0.832 | 0.637 | 0.645 | 1 |
| 10 | Idelalisib | 0.555 | 0.746 | 0.759 | 0.976 | 0.254 | 0.698 | 1 |

## clinical-trial agent  (36 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | OTX-015 | 0.719 | 0.880 | 0.779 | 0.998 | 0.920 | 0.541 | 0 |
| 2 | AZD1480 | 0.659 | 0.686 | 0.680 | 0.995 | 0.973 | 0.635 | 0 |
| 3 | Saracatinib (AZD0530) | 0.657 | 0.775 | 0.872 | 0.921 | 0.783 | 0.642 | 0 |
| 4 | Cediranib (AZD2171) | 0.569 | 0.793 | 0.748 | 0.982 | 0.563 | 0.434 | 0 |
| 5 | AT7519 | 0.563 | 0.819 | 0.787 | 0.996 | 0.893 | 0.328 | 2 |
| 6 | GDC-0941 | 0.562 | 0.893 | 0.825 | 0.977 | 0.254 | 0.569 | 2 |
| 7 | Elesclomol | 0.547 | 0.737 | 0.800 | 0.805 | 0.705 | 0.388 | 1 |
| 8 | 17-AAG (Tanespimycin) | 0.532 | 0.940 | 0.894 | 0.976 | 0.254 | 0.233 | 1 |
| 9 | CI-1040 (PD184352) | 0.522 | 0.851 | 0.819 | 0.971 | 0.254 | 0.460 | 0 |
| 10 | INK-128 | 0.503 | 0.781 | 0.793 | 0.977 | 0.254 | 0.518 | 1 |

## research-only compound  (40 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | PRT062607 | 0.703 | 0.842 | 0.804 | 0.991 | 0.773 | 0.853 | 0 |
| 2 | JAK Inhibitor I | 0.595 | 0.716 | 0.700 | 0.995 | 0.815 | 0.491 | 0 |
| 3 | JQ1 | 0.565 | 0.819 | 0.763 | 0.984 | 0.484 | 0.541 | 0 |
| 4 | AGI-5198 | 0.535 | 0.546 | 0.649 | 0.632 | 1.000 | 0.683 | 1 |
| 5 | PP242 | 0.528 | 0.878 | 0.816 | 0.977 | 0.254 | 0.518 | 1 |
| 6 | Doramapimod (BIRB 796) | 0.496 | 0.888 | 0.841 | 0.977 | 0.271 | 0.457 | 2 |
| 7 | Nutlin 3a | 0.441 | 0.417 | 0.786 | 0.493 | 1.000 | 0.433 | 3 |
| 8 | TG100-115 | 0.438 | 0.650 | 0.689 | 0.977 | 0.254 | 0.635 | 1 |
| 9 | JNJ-28312141 | 0.410 | 0.574 | 0.842 | 0.683 | 0.516 | 0.630 | 1 |
| 10 | SU11274 | 0.391 | 0.460 | 0.809 | 0.365 | 0.729 | 0.426 | 0 |

## Combination hypotheses (complementary cell-state coverage)

*BeatAML2 measured single agents only; no synergy is claimed.*

| pair | pathways | coverage A | coverage B | union | gain |
|---|---|---|---|---|---|
| Crizotinib (PF-2341066) + Midostaurin | RTK / FLT3 | 0.47 | 0.24 | 0.69 | +0.22 |
| Cytarabine + Lenalidomide | chemotherapy / epigenetic | 0.69 | 0.78 | 1.00 | +0.22 |
| Lenalidomide + Midostaurin | epigenetic / FLT3 | 0.78 | 0.24 | 1.00 | +0.22 |
| Midostaurin + Palbociclib | FLT3 / cell_cycle | 0.24 | 0.47 | 0.69 | +0.22 |
| Azacytidine + Midostaurin | epigenetic / FLT3 | 0.17 | 0.24 | 0.40 | +0.16 |
| Gefitinib + Midostaurin | RTK / FLT3 | 0.17 | 0.24 | 0.40 | +0.16 |
