# IRS: Iterative Reference Selection for Microbiome Data

IRS selects stable reference taxa and estimates reusable sample normalization
factors for microbiome differential abundance analysis. Fit IRS once, inspect
the selected references, and continue with either the built-in Wilcoxon test or
a supported downstream DAA package.

Counts are supplied as a taxa-by-samples matrix with unique taxon and sample
names.

## Installation

```r
# install.packages("remotes")
remotes::install_github("yimshi/IRS", build_vignettes = TRUE)
```

## Basic workflow

```r
library(IRS)

fit <- fit_irs(
  counts = counts,
  meta_dat = metadata,
  predictor = "group",
  seed = 123
)

fit
reference_taxa(fit)
diagnostics(fit)

irs_counts <- normalized_counts(fit, counts)
result <- irs_wilcox(
  fit,
  counts,
  group = "group",
  metadata = metadata
)
```

For sample *i*, IRS uses the final reference total
`T_i^R = sum(counts[reference_taxa, i])`. The stored size factors are centered
to a geometric mean of one, which keeps normalized values on a convenient
count-like scale.

## Downstream integrations

| Method | Helper | Integration |
|---|---|---|
| DESeq2 | `inject_deseq2()` | IRS size factors |
| edgeR | `inject_edger()` | IRS effective library sizes |
| MaAsLin2 LM | `prepare_maaslin2()` | IRS-normalized input |
| ZicoSeq | `zicoseq_irs()` | IRS references and totals |
| LDM frequency scale | `prepare_ldm()` | IRS-normalized input |
| Wilcoxon | `irs_wilcox()` | complete two-group workflow |
| DACOMP | `reference_for_dacomp()` | IRS-selected references |
| ALDEx2 | `denom_for_aldex2()` | IRS-selected denominator taxa |

For example, IRS factors can be inserted into DESeq2 without changing the
integer count matrix:

```r
dds <- DESeq2::DESeqDataSetFromMatrix(counts, metadata, ~ group)
dds <- inject_deseq2(dds, fit)
dds <- DESeq2::DESeq(dds)
result <- DESeq2::results(dds, contrast = c("group", "case", "control"))
```

## Vignettes and example data

```r
vignette("downstream-daa", package = "IRS")
vignette("vandeputte2017", package = "IRS")

data("vandeputte2017", package = "IRS")
dim(vandeputte2017$counts)
table(vandeputte2017$sample_data$group)
```

The first vignette is a runnable introduction and includes examples for every
supported adapter. The second uses the Vandeputte 2017 Crohn's disease cohort
and compares IRS-normalized results with flow-cytometry-based quantitative
microbiome profiles.

The original `select_reference_irs()` interface remains available for code
written against earlier IRS versions.
