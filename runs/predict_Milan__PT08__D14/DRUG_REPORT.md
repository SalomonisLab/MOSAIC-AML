# MOSAIC-Rx drug prioritisation

*Predicted **ex-vivo** sensitivity from the BeatAML2 functional screen. This is a prioritisation for trial matching or laboratory validation, **not** a treatment recommendation.*

- inhibitors modelled: **118**, reported 112, abstained 6
- genes shared with the model: 14237
- distance from the BeatAML training distribution: 98th percentile
- differentiation axis (primitive - mature): 63th percentile of BeatAML
- cell states scored: 2

## approved in AML  (6 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Venetoclax | 0.711 | 0.987 | 0.966 | — | — | 0.823 | 1 |
| 2 | Midostaurin | 0.397 | 0.635 | 0.803 | — | — | 0.379 | 1 |
| 3 | Azacytidine | 0.241 | 0.380 | 0.646 | — | — | 0.485 | 2 |
| 4 | Cytarabine | 0.181 | 0.395 | 0.824 | — | — | 0.099 | 2 |
| 5 | Quizartinib (AC220) | 0.107 | 0.070 | 0.855 | — | — | 0.197 | 1 |
| 6 | Gilteritinib | 0.091 | 0.048 | 0.850 | — | — | 0.187 | 1 |

## approved, other indication (off-label in AML)  (35 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Neratinib (HKI-272) | 0.492 | 0.814 | 0.716 | — | — | 0.673 | 1 |
| 2 | Erlotinib | 0.406 | 0.776 | 0.834 | — | — | 0.357 | 2 |
| 3 | Lapatinib | 0.395 | 0.650 | 0.734 | — | — | 0.673 | 1 |
| 4 | Trametinib (GSK1120212) | 0.393 | 0.664 | 0.903 | — | — | 0.539 | 1 |
| 5 | Vandetanib (ZD6474) | 0.381 | 0.757 | 0.825 | — | — | 0.309 | 2 |
| 6 | Afatinib (BIBW-2992) | 0.377 | 0.624 | 0.718 | — | — | 0.673 | 1 |
| 7 | Ibrutinib (PCI-32765) | 0.373 | 0.754 | 0.825 | — | — | 0.287 | 2 |
| 8 | Axitinib (AG-013736) | 0.371 | 0.719 | 0.862 | — | — | 0.341 | 1 |
| 9 | Nilotinib | 0.358 | 0.739 | 0.690 | — | — | 0.334 | 2 |
| 10 | Rapamycin | 0.338 | 0.553 | 0.922 | — | — | 0.579 | 1 |

## clinical-trial agent  (36 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | LY-333531 | 0.466 | 0.886 | 0.871 | — | — | 0.567 | 1 |
| 2 | Canertinib (CI-1033) | 0.461 | 0.855 | 0.773 | — | — | 0.673 | 1 |
| 3 | RAF265 (CHIR-265) | 0.406 | 0.848 | 0.831 | — | — | 0.437 | 2 |
| 4 | CI-1040 (PD184352) | 0.302 | 0.637 | 0.819 | — | — | 0.539 | 1 |
| 5 | AT7519 | 0.284 | 0.771 | 0.787 | — | — | 0.159 | 2 |
| 6 | Masitinib (AB-1010) | 0.275 | 0.662 | 0.684 | — | — | 0.441 | 1 |
| 7 | MGCD-265 | 0.247 | 0.652 | 0.843 | — | — | 0.269 | 1 |
| 8 | VX-745 | 0.244 | 0.637 | 0.666 | — | — | 0.388 | 1 |
| 9 | BEZ235 | 0.242 | 0.558 | 0.765 | — | — | 0.517 | 1 |
| 10 | Alisertib (MLN8237) | 0.237 | 0.764 | 0.728 | — | — | 0.019 | 2 |

## research-only compound  (35 considered)

| # | inhibitor | utility | P(sensitive) | model AUROC | blast cov. | LSC cov. | mech | challenges |
|---|---|---|---|---|---|---|---|---|
| 1 | Nutlin 3a | 0.462 | 0.796 | 0.786 | — | — | 0.920 | 2 |
| 2 | ABT-737 | 0.456 | 0.832 | 0.763 | — | — | 0.833 | 1 |
| 3 | NF-kB Activation Inhibitor | 0.407 | 0.916 | 0.865 | — | — | 0.558 | 2 |
| 4 | JQ1 | 0.272 | 0.717 | 0.763 | — | — | 0.365 | 1 |
| 5 | PHA-665752 | 0.271 | 0.757 | 0.798 | — | — | 0.236 | 2 |
| 6 | KU-55933 | 0.253 | 0.707 | 0.737 | — | — | 0.516 | 2 |
| 7 | NVP-ADW742 | 0.247 | 0.694 | 0.723 | — | — | 0.531 | 2 |
| 8 | GSK-1838705A | 0.246 | 0.761 | 0.744 | — | — | 0.355 | 3 |
| 9 | PD173955 | 0.221 | 0.578 | 0.794 | — | — | 0.470 | 1 |
| 10 | SR9011 | 0.208 | 0.647 | 0.721 | — | — | 0.508 | 2 |
