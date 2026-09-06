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
})

test_that("edgeR effective library sizes are proportional to IRS totals", {
  skip_if_not_installed("edgeR")
  data <- make_irs_test_data()
  fit <- make_irs_test_fit(data)
  y <- edgeR::DGEList(data$counts)
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

  ldm <- prepare_ldm(fit, data$counts)
  expect_false(ldm$scale.otu.table)
  expect_true(ldm$freq.scale.only)

  dacomp_ref <- reference_for_dacomp(fit, t(data$counts))
  expect_identical(
    colnames(t(data$counts))[unname(dacomp_ref)], reference_taxa(fit)
  )

  aldex_ref <- denom_for_aldex2(fit, data$counts)
  expect_identical(rownames(data$counts)[unname(aldex_ref)], reference_taxa(fit))
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
})
