# IRS 0.2.0

* Added `fit_irs()` and the `irs_fit` class with reference totals, centered size
  factors, log offsets, convergence traces, sample diagnostics, and optional
  count storage.
* Added accessors for reference taxa, indices, totals, size factors, offsets,
  normalized counts, and diagnostics.
* Added native normalization injection for DESeq2 and edgeR.
* Added prepared-input helpers for the MaAsLin2 linear model and LDM
  frequency-scale analysis.
* Added IRS reference-set helpers for DACOMP and ALDEx2.
* Added taxon-wise IRS-normalized Wilcoxon testing.
* Added `zicoseq_irs()`, derived from GUniFrac 1.9 ZicoSeq under GPL-3, with a
  fixed named IRS reference set and original-count reference denominators.
* Changed the package license from MIT to GPL-3 to support distribution of the
  ZicoSeq-derived implementation.
* Retained `select_reference_irs()` for backward compatibility and added stricter
  validation and convergence reporting through `fit_irs()`.
* Added a runnable getting-started vignette, a Vandeputte 2017 case-study
  vignette, and precomputed summaries of paired downstream DAA benchmarks.
* Added the `vandeputte2017` example data with genus-level 16S counts,
  Crohn's-disease status, and flow-cytometry microbial-load measurements.
