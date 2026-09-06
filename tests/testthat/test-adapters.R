test_that("DESeq2 receives IRS size factors without changing counts", {
  skip_if_not_installed("DESeq2")
  data <- make_irs_test_data()
  fit <- make_irs_test_fit(data)
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = data$counts,
    colData = data$metadata,
    design = ~ group
  )
  original <- DESeq2::counts(dds)
  dds <- inject_deseq2(dds, fit)

  expect_equal(DESeq2::counts(dds), original)
  expect_equal(
    unname(DESeq2::sizeFactors(dds)),
    unname(size_factors(fit)),
    tolerance = 1e-12
  )

  injected_size_factors <- DESeq2::sizeFactors(dds)
  dds <- DESeq2::DESeq(dds, fitType = "mean", quiet = TRUE)
  expect_equal(
    unname(DESeq2::sizeFactors(dds)),
    unname(injected_size_factors),
    tolerance = 1e-12
  )
})

test_that("edgeR effective library sizes are proportional to IRS totals", {
  skip_if_not_installed("edgeR")
  data <- make_irs_test_data()
  fit <- make_irs_test_fit(data)
  y <- edgeR::DGEList(data$counts)
  design <- stats::model.matrix(~ group, data$metadata)
  keep <- edgeR::filterByExpr(y, design = design)
  y <- y[keep, , keep.lib.sizes = FALSE]
  original <- y$counts
  y <- inject_edger(y, fit)
  effective <- y$samples$lib.size * y$samples$norm.factors

  expect_equal(y$counts, original)
  expect_lt(stats::sd(log(effective / reference_totals(fit))), 1e-12)
  expect_equal(exp(mean(log(y$samples$norm.factors))), 1, tolerance = 1e-12)
})

test_that("matrix and reference-set adapters preserve orientation and names", {
  data <- make_irs_test_data()
  fit <- make_irs_test_fit(data)

  m2 <- prepare_maaslin2(fit, data$counts)
  expect_identical(rownames(m2$input_data), colnames(data$counts))
  expect_identical(colnames(m2$input_data), rownames(data$counts))
  expect_identical(m2$normalization, "NONE")
  expect_identical(m2$analysis_method, "LM")

  ldm <- prepare_ldm(fit, data$counts)
  expect_false(ldm$scale.otu.table)
  expect_true(ldm$freq.scale.only)
  expect_false(ldm$comp.anal)
  expect_false(all(abs(rowSums(ldm$otu_table) - 1) < 1e-12))

  dacomp_ref <- reference_for_dacomp(fit, t(data$counts))
  expect_identical(
    colnames(t(data$counts))[unname(dacomp_ref)], reference_taxa(fit)
  )

  aldex_ref <- denom_for_aldex2(fit, data$counts)
  expect_identical(rownames(data$counts)[unname(aldex_ref)], reference_taxa(fit))

  reads_with_zero <- rbind(all_zero = 0, data$counts)
  shifted_aldex_ref <- denom_for_aldex2(fit, reads_with_zero)
  kept_names <- rownames(reads_with_zero)[rowSums(reads_with_zero) > 0]
  expect_identical(
    kept_names[unname(shifted_aldex_ref)], reference_taxa(fit)
  )
})

test_that("IRS Wilcoxon returns named adjusted results", {
  data <- make_irs_test_data()
  fit <- make_irs_test_fit(data)
  result <- irs_wilcox(
    fit, data$counts, group = "group", metadata = data$metadata
  )

  expect_identical(result$taxon, rownames(data$counts))
  expect_true(all(result$p_value >= 0 & result$p_value <= 1))
  expect_true(all(result$p_adjusted >= 0 & result$p_adjusted <= 1))

  three_groups <- rep(c("a", "b", "c"), length.out = ncol(data$counts))
  expect_error(
    irs_wilcox(fit, data$counts, group = three_groups),
    "exactly two groups"
  )
})
