# MOSAIC-Rx drug prioritisation

*Predicted **ex-vivo** sensitivity from the BeatAML2 functional screen. This is a prioritisation for trial matching or laboratory validation, **not** a treatment recommendation.*

- inhibitors modelled: **118**, reported 118, abstained 0
- genes shared with the model: 14237
- distance from the BeatAML training distribution: 17th percentile
- differentiation axis (primitive - mature): 36th percentile of BeatAML
- cell states scored: 44

## approved in AML  (6 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Quizartinib (AC220) | 0.851 | 0.877 | 0.855 | 1.000 | 1.000 | 0.572 | 0 |
| 2 | Midostaurin | 0.744 | 0.792 | 0.803 | 0.881 | 0.854 | 0.540 | 0 |
| 3 | Gilteritinib | 0.659 | 0.710 | 0.850 | 0.653 | 0.739 | 0.621 | 1 |
| 4 | Cytarabine | 0.341 | 0.487 | 0.824 | 0.207 | 0.000 | 0.452 | 1 |
| 5 | Venetoclax | 0.289 | 0.047 | 0.966 | 0.169 | 0.169 | 0.625 | 2 |
| 6 | Azacytidine | 0.223 | 0.252 | 0.646 | 0.000 | 0.000 | 0.450 | 2 |

## approved, other indication (off-label in AML)  (36 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Dasatinib | 0.870 | 0.960 | 0.931 | 0.996 | 1.000 | 0.739 | 0 |
| 2 | Axitinib (AG-013736) | 0.845 | 0.903 | 0.862 | 1.000 | 1.000 | 0.773 | 0 |
| 3 | Trametinib (GSK1120212) | 0.843 | 0.928 | 0.903 | 0.981 | 1.000 | 0.724 | 0 |
| 4 | Vargetef | 0.838 | 0.882 | 0.827 | 1.000 | 1.000 | 0.803 | 0 |
| 5 | Lenvatinib | 0.835 | 0.874 | 0.861 | 0.992 | 1.000 | 0.803 | 0 |
| 6 | Tivozanib (AV-951) | 0.821 | 0.951 | 0.941 | 0.887 | 0.854 | 0.801 | 0 |
| 7 | Bosutinib (SKI-606) | 0.806 | 0.906 | 0.862 | 0.910 | 0.854 | 0.862 | 0 |
| 8 | Regorafenib (BAY 73-4506) | 0.797 | 0.836 | 0.834 | 0.945 | 0.977 | 0.763 | 0 |
| 9 | Pazopanib (GW786034) | 0.792 | 0.891 | 0.885 | 0.929 | 0.854 | 0.770 | 0 |
| 10 | Panobinostat | 0.790 | 0.867 | 0.860 | 0.996 | 0.992 | 0.580 | 0 |

## clinical-trial agent  (36 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Saracatinib (AZD0530) | 0.830 | 0.931 | 0.872 | 1.000 | 1.000 | 0.861 | 0 |
| 2 | Linifanib (ABT-869) | 0.805 | 0.930 | 0.903 | 0.995 | 1.000 | 0.702 | 0 |
| 3 | Entospletinib (GS-9973) | 0.788 | 0.943 | 0.872 | 0.985 | 1.000 | 0.628 | 0 |
| 4 | CI-1040 (PD184352) | 0.764 | 0.857 | 0.819 | 0.981 | 1.000 | 0.724 | 0 |
| 5 | Crenolanib | 0.763 | 0.902 | 0.879 | 0.970 | 1.000 | 0.562 | 0 |
| 6 | YM-155 | 0.732 | 0.816 | 0.754 | 1.000 | 1.000 | 0.635 | 0 |
| 7 | Cediranib (AZD2171) | 0.718 | 0.761 | 0.748 | 0.985 | 0.969 | 0.773 | 0 |
| 8 | LY-333531 | 0.687 | 0.880 | 0.871 | 0.857 | 0.736 | 0.707 | 0 |
| 9 | Motesanib (AMG-706) | 0.687 | 0.776 | 0.824 | 0.868 | 0.854 | 0.791 | 0 |
| 10 | OTX-015 | 0.675 | 0.762 | 0.779 | 0.944 | 1.000 | 0.474 | 0 |

## research-only compound  (40 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | PRT062607 | 0.729 | 0.849 | 0.804 | 0.999 | 1.000 | 0.628 | 0 |
| 2 | JNJ-28312141 | 0.718 | 0.861 | 0.842 | 0.997 | 1.000 | 0.689 | 1 |
| 3 | KI20227 | 0.659 | 0.844 | 0.808 | 0.929 | 0.854 | 0.696 | 1 |
| 4 | PD173955 | 0.634 | 0.733 | 0.794 | 0.868 | 0.739 | 0.862 | 0 |
| 5 | Doramapimod (BIRB 796) | 0.615 | 0.859 | 0.841 | 0.863 | 0.708 | 0.684 | 1 |
| 6 | JNJ-7706621 | 0.590 | 0.672 | 0.663 | 1.000 | 1.000 | 0.305 | 0 |
| 7 | PHA-665752 | 0.528 | 0.581 | 0.798 | 0.557 | 0.854 | 0.744 | 0 |
| 8 | Nutlin 3a | 0.491 | 0.472 | 0.786 | 0.575 | 0.900 | 0.640 | 2 |
| 9 | NF-kB Activation Inhibitor | 0.490 | 0.461 | 0.865 | 0.586 | 0.969 | 0.711 | 2 |
| 10 | GSK-1838705A | 0.489 | 0.614 | 0.744 | 0.596 | 0.854 | 0.634 | 2 |

## Combination hypotheses (complementary cell-state coverage)

*BeatAML2 measured single agents only; no synergy is claimed.*

| pair | pathways | coverage A | coverage B | union | gain |
|---|---|---|---|---|---|
| Bortezomib (Velcade) + Lenalidomide | proteostasis / epigenetic | 0.49 | 0.57 | 0.93 | +0.35 |
| Ibrutinib (PCI-32765) + Lenalidomide | immune_signalling / epigenetic | 0.57 | 0.57 | 0.93 | +0.35 |
| Idelalisib + Lenalidomide | PI3K_AKT_mTOR / epigenetic | 0.57 | 0.57 | 0.93 | +0.35 |
| Crizotinib (PF-2341066) + Lenalidomide | RTK / epigenetic | 0.60 | 0.57 | 0.93 | +0.33 |
| Bortezomib (Velcade) + Erlotinib | proteostasis / RTK | 0.49 | 0.58 | 0.88 | +0.30 |
| Erlotinib + Ibrutinib (PCI-32765) | RTK / immune_signalling | 0.58 | 0.57 | 0.88 | +0.30 |
| Erlotinib + Idelalisib | RTK / PI3K_AKT_mTOR | 0.58 | 0.57 | 0.88 | +0.30 |
| Bortezomib (Velcade) + Palbociclib | proteostasis / cell_cycle | 0.49 | 0.74 | 1.00 | +0.26 |
