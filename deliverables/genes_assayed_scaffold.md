# "Genes assayed" coverage matrix — scaffold

**The problem this fixes.** In the single-cell cohorts, a gene a center never sequenced is recorded as
`absent`. So an `absent` label conflates two very different things:

- *assayed and wild-type* — a true negative, and
- *not assayed* — **no information**, wrongly counted as a true negative.

When the model calls such a sample `present`, it is scored as a **false positive** even though the
"truth" was never established. This **understates specificity and precision**. The fix is a per-sample
record of which genes were actually tested.

## The format (one row per sample, one column per gene/category)

`labels/genes_assayed.tsv` — values: `1` assayed, `0` not assayed (blank = unknown).

```
sample_key                assay_type            TP53  FLT3  NPM1  DNMT3A  IDH1  IDH2  inv16  del5  ...
CCHMC::1009_AfInv16_29M    karyotype+panel-27    1     1     1     1       1     1     1      1     ...
Colorado::AML-04           WES                   1     1     1     1       1     1     1      1     ...
NYU-1::AML0024             targeted-54           1     1     1     1       1     1     0      0     ...
...
```

- **Karyotype** covers the cytogenetic events (inv16, del5, del7, trisomy8, complex, KMT2A) — reliably assayed for essentially every sample.
- **Targeted panels** cover a defined gene list (varies by center); genes outside the panel → `0`.
- **WES** covers all SNV genes → `1` across the board.

## How it changes the metrics (the correction)

For each mutation, restrict the negatives to **assayed-wild-type only**:

```
specificity* = TN_assayed / (TN_assayed + FP_assayed)      # drop not-assayed samples from the negative set
precision*    = TP / (TP + FP_assayed)
```

Samples where the gene was **not assayed** are excluded from that mutation's negatives (they carry no
truth), exactly as we already do for the MDS empty-`Group` samples. Expected effect: specificity and
precision **rise** (some current "false positives" were never really negatives), while sensitivity is
unchanged (positives are still positives).

## Status

- ◑ **Cytogenetics**: assayable now — karyotype is ~universal, so the 6 cytogenetic columns are `1` for all sc samples with a reported karyotype.
- ☐ **SNV panels**: needs the **per-center panel gene lists** (and which samples got WES vs targeted). This is a metadata request to the 6 contributing centers (CCHMC, Colorado, Columbia, Milan, NYU, WashU), not a computation.

Until the panel metadata arrives, single-cell **specificity is a lower bound** — the true value is at
least as good as reported. That statement itself is worth putting in the manuscript.
