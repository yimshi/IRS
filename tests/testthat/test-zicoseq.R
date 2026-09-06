test_that("zicoseq_irs uses fixed named IRS references", {
  skip_if_not_installed("vegan")
  data <- make_irs_test_data(n_taxa = 16L, n_samples = 50L)
  fit <- make_irs_test_fit(data)

  set.seed(99)
  result <- zicoseq_irs(
    meta.dat = data$metadata,
    feature.dat = data$counts,
    grp.name = "group",
    is.winsor = FALSE,
    is.post.sample = FALSE,
    perm.no = 99,
    verbose = FALSE,
    irs_fit = fit
  )

  expect_identical(result$ref.features, reference_taxa(fit))
  expect_identical(result$irs$denominator_source, "original_counts")
  expect_equal(
    unname(result$irs$raw_reference_total),
    unname(colSums(data$counts[reference_taxa(fit), , drop = FALSE]))
  )
  expect_named(result$p.raw, rownames(result$F0))
  expect_length(result$p.adj.fdr, nrow(data$counts))
})

test_that("zicoseq_irs reference totals do not change with feature filtering", {
  skip_if_not_installed("vegan")
  data <- make_irs_test_data(n_taxa = 16L, n_samples = 50L)
  fit <- make_irs_test_fit(data)
  expected <- colSums(data$counts[reference_taxa(fit), , drop = FALSE])

  set.seed(100)
  result <- zicoseq_irs(
    data$metadata,
    data$counts,
    "group",
    mean.abund.filter = 0.03,
    is.post.sample = FALSE,
    perm.no = 99,
    verbose = FALSE,
    irs_fit = fit
  )
  expect_equal(unname(result$irs$raw_reference_total), unname(expected))
})

test_that("zicoseq_irs retains posterior sampling and winsorization paths", {
  skip_if_not_installed("vegan")
  skip_if_not_installed("rmutil")
  data <- make_irs_test_data(n_taxa = 10L, n_samples = 50L, seed = 91L)
  fit <- make_irs_test_fit(data)

  set.seed(101)
  result <- zicoseq_irs(
    data$metadata,
    data$counts,
    "group",
    is.winsor = TRUE,
    outlier.pct = 0.01,
    is.post.sample = TRUE,
    post.sample.no = 2,
    perm.no = 99,
    verbose = FALSE,
    irs_fit = fit
  )
  expect_length(result$p.raw, nrow(data$counts))
  expect_equal(
    unname(result$irs$raw_reference_total),
    unname(reference_totals(fit))
  )
})
