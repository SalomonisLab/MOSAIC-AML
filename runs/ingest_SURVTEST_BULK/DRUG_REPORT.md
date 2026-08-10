# COMPASS-AML drug prioritisation

*Predicted **ex-vivo** sensitivity from the BeatAML2 functional screen. This is a prioritisation for trial matching or laboratory validation, **not** a treatment recommendation.*

- inhibitors modelled: **118**, reported 113, abstained 5
- genes shared with the model: 14163
- distance from the BeatAML training distribution: 16th percentile
- differentiation axis (primitive - mature): 80th percentile of BeatAML
- cell states scored: bulk input, none

## approved in AML  (6 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Cytarabine | 0.755 | 0.729 | 0.824 | — | — | 0.822 | 1 |
| 2 | Venetoclax | 0.652 | 0.759 | 0.966 | — | — | 0.438 | 0 |
| 3 | Azacytidine | 0.523 | 0.262 | 0.646 | — | — | 0.915 | 2 |
| 4 | Gilteritinib | 0.432 | 0.373 | 0.850 | — | — | 0.635 | 0 |
| 5 | Midostaurin | 0.402 | 0.345 | 0.803 | — | — | 0.569 | 0 |
| 6 | Quizartinib (AC220) | 0.341 | 0.151 | 0.855 | — | — | 0.606 | 1 |

## approved, other indication (off-label in AML)  (35 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Pazopanib (GW786034) | 0.663 | 0.739 | 0.885 | — | — | 0.687 | 0 |
| 2 | Bosutinib (SKI-606) | 0.657 | 0.736 | 0.862 | — | — | 0.682 | 0 |
| 3 | Afatinib (BIBW-2992) | 0.597 | 0.665 | 0.718 | — | — | 0.691 | 0 |
| 4 | Selumetinib (AZD6244) | 0.561 | 0.597 | 0.886 | — | — | 0.622 | 0 |
| 5 | Erlotinib | 0.547 | 0.543 | 0.834 | — | — | 0.722 | 0 |
| 6 | Sorafenib | 0.539 | 0.716 | 0.919 | — | — | 0.586 | 0 |
| 7 | Neratinib (HKI-272) | 0.524 | 0.546 | 0.716 | — | — | 0.691 | 0 |
| 8 | Palbociclib | 0.515 | 0.364 | 0.782 | — | — | 0.957 | 1 |
| 9 | Vandetanib (ZD6474) | 0.513 | 0.502 | 0.825 | — | — | 0.693 | 0 |
| 10 | Trametinib (GSK1120212) | 0.505 | 0.499 | 0.903 | — | — | 0.622 | 0 |

## clinical-trial agent  (35 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | YM-155 | 0.566 | 0.731 | 0.754 | — | — | 0.639 | 0 |
| 2 | KW-2449 | 0.533 | 0.601 | 0.913 | — | — | 0.733 | 0 |
| 3 | Cediranib (AZD2171) | 0.532 | 0.653 | 0.748 | — | — | 0.692 | 0 |
| 4 | OTX-015 | 0.506 | 0.552 | 0.779 | — | — | 0.831 | 0 |
| 5 | CI-1040 (PD184352) | 0.499 | 0.613 | 0.819 | — | — | 0.622 | 0 |
| 6 | Pelitinib (EKB-569) | 0.459 | 0.484 | 0.893 | — | — | 0.722 | 0 |
| 7 | RAF265 (CHIR-265) | 0.455 | 0.518 | 0.831 | — | — | 0.669 | 0 |
| 8 | Alisertib (MLN8237) | 0.448 | 0.634 | 0.728 | — | — | 0.774 | 0 |
| 9 | Dovitinib (CHIR-258) | 0.444 | 0.441 | 0.917 | — | — | 0.715 | 0 |
| 10 | GSK690693 | 0.435 | 0.553 | 0.658 | — | — | 0.600 | 0 |

## research-only compound  (37 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | NVP-TAE684 | 0.642 | 0.829 | 0.907 | — | — | 0.722 | 0 |
| 2 | GSK-1838705A | 0.519 | 0.739 | 0.744 | — | — | 0.739 | 1 |
| 3 | PHA-665752 | 0.464 | 0.599 | 0.798 | — | — | 0.624 | 0 |
| 4 | GW-2580 | 0.461 | 0.567 | 0.863 | — | — | 0.655 | 0 |
| 5 | NVP-ADW742 | 0.456 | 0.633 | 0.723 | — | — | 0.754 | 1 |
| 6 | KU-55933 | 0.433 | 0.561 | 0.737 | — | — | 0.827 | 1 |
| 7 | PP242 | 0.428 | 0.505 | 0.816 | — | — | 0.698 | 0 |
| 8 | SR9011 | 0.408 | 0.581 | 0.721 | — | — | 0.701 | 1 |
| 9 | GSK-1904529A | 0.406 | 0.562 | 0.696 | — | — | 0.739 | 1 |
| 10 | PD173955 | 0.390 | 0.430 | 0.794 | — | — | 0.682 | 0 |
