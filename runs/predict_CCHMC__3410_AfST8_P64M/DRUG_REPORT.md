# MOSAIC-Rx drug prioritisation

*Predicted **ex-vivo** sensitivity from the BeatAML2 functional screen. This is a prioritisation for trial matching or laboratory validation, **not** a treatment recommendation.*

- inhibitors modelled: **118**, reported 118, abstained 0
- genes shared with the model: 14237
- distance from the BeatAML training distribution: 70th percentile
- differentiation axis (primitive - mature): 80th percentile of BeatAML
- cell states scored: 34

## approved in AML  (6 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Gilteritinib | 0.766 | 0.912 | 0.850 | 0.984 | 1.000 | 0.508 | 1 |
| 2 | Quizartinib (AC220) | 0.744 | 0.857 | 0.855 | 0.982 | 1.000 | 0.505 | 0 |
| 3 | Midostaurin | 0.694 | 0.790 | 0.803 | 0.994 | 1.000 | 0.390 | 1 |
| 4 | Venetoclax | 0.692 | 0.889 | 0.966 | 0.986 | 0.992 | 0.055 | 1 |
| 5 | Azacytidine | 0.148 | 0.372 | 0.646 | 0.001 | 0.001 | 0.130 | 2 |
| 6 | Cytarabine | 0.139 | 0.120 | 0.824 | 0.001 | 0.000 | 0.336 | 2 |

## approved, other indication (off-label in AML)  (36 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Ponatinib (AP24534) | 0.706 | 0.902 | 0.952 | 0.940 | 0.953 | 0.522 | 0 |
| 2 | Gefitinib | 0.667 | 0.809 | 0.794 | 0.994 | 1.000 | 0.486 | 0 |
| 3 | Regorafenib (BAY 73-4506) | 0.637 | 0.743 | 0.834 | 0.972 | 0.998 | 0.479 | 0 |
| 4 | Palbociclib | 0.628 | 0.772 | 0.782 | 0.993 | 0.999 | 0.369 | 1 |
| 5 | Lenvatinib | 0.623 | 0.693 | 0.861 | 0.968 | 1.000 | 0.513 | 0 |
| 6 | Lapatinib | 0.621 | 0.635 | 0.734 | 0.981 | 0.999 | 0.685 | 0 |
| 7 | Vargetef | 0.616 | 0.676 | 0.827 | 0.984 | 1.000 | 0.500 | 0 |
| 8 | Neratinib (HKI-272) | 0.613 | 0.629 | 0.716 | 0.981 | 1.000 | 0.685 | 0 |
| 9 | CYT387 | 0.603 | 0.701 | 0.703 | 0.997 | 0.999 | 0.430 | 0 |
| 10 | Ruxolitinib (INCB018424) | 0.550 | 0.575 | 0.713 | 0.993 | 0.999 | 0.419 | 0 |

## clinical-trial agent  (36 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Entospletinib (GS-9973) | 0.677 | 0.905 | 0.872 | 0.983 | 1.000 | 0.524 | 0 |
| 2 | Pelitinib (EKB-569) | 0.622 | 0.782 | 0.893 | 0.980 | 1.000 | 0.486 | 0 |
| 3 | GDC-0941 | 0.617 | 0.797 | 0.825 | 0.986 | 1.000 | 0.475 | 1 |
| 4 | LY-333531 | 0.614 | 0.880 | 0.871 | 0.993 | 0.999 | 0.195 | 1 |
| 5 | Canertinib (CI-1033) | 0.606 | 0.698 | 0.773 | 0.980 | 1.000 | 0.685 | 0 |
| 6 | Volasertib (BI-6727) | 0.561 | 0.707 | 0.668 | 0.993 | 1.000 | 0.466 | 0 |
| 7 | Dovitinib (CHIR-258) | 0.557 | 0.658 | 0.917 | 0.898 | 1.000 | 0.477 | 0 |
| 8 | Flavopiridol | 0.550 | 0.714 | 0.695 | 0.997 | 1.000 | 0.342 | 0 |
| 9 | AT7519 | 0.550 | 0.694 | 0.787 | 0.988 | 0.999 | 0.336 | 1 |
| 10 | 17-AAG (Tanespimycin) | 0.537 | 0.565 | 0.894 | 0.920 | 0.951 | 0.661 | 0 |

## research-only compound  (40 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | NF-kB Activation Inhibitor | 0.651 | 0.934 | 0.865 | 0.985 | 1.000 | 0.601 | 1 |
| 2 | PRT062607 | 0.580 | 0.724 | 0.804 | 0.988 | 0.999 | 0.524 | 0 |
| 3 | JAK Inhibitor I | 0.553 | 0.706 | 0.700 | 0.996 | 1.000 | 0.473 | 0 |
| 4 | PI-103 | 0.547 | 0.709 | 0.757 | 0.945 | 0.953 | 0.529 | 0 |
| 5 | PP242 | 0.529 | 0.657 | 0.816 | 0.945 | 0.953 | 0.510 | 0 |
| 6 | NVP-TAE684 | 0.527 | 0.643 | 0.907 | 0.881 | 0.979 | 0.486 | 0 |
| 7 | Bay 11-7085 | 0.522 | 0.667 | 0.705 | 0.993 | 1.000 | 0.593 | 1 |
| 8 | NVP-ADW742 | 0.519 | 0.768 | 0.723 | 0.987 | 1.000 | 0.322 | 2 |
| 9 | ABT-737 | 0.513 | 0.744 | 0.763 | 0.978 | 0.998 | 0.142 | 2 |
| 10 | SU11274 | 0.491 | 0.741 | 0.809 | 0.973 | 0.997 | 0.522 | 0 |
