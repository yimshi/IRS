test_that("fit_irs returns centered factors and exact IRS ratios", {
  data <- make_irs_test_data()
  fit <- make_irs_test_fit(data)

  expect_s3_class(fit, "irs_fit")
  expect_true(fit$converged)
  expect_equal(exp(mean(log(size_factors(fit)))), 1, tolerance = 1e-12)
  expect_identical(reference_taxa(fit), rownames(data$counts)[fit$reference])

  normalized <- normalized_counts(fit, data$counts)
  raw_ratio <- sweep(data$counts, 2, reference_totals(fit), "/")
  expected <- raw_ratio * exp(mean(log(reference_totals(fit))))
  expect_equal(unclass(normalized), unclass(expected), tolerance = 1e-12)
})

test_that("fit_irs aligns metadata and exposes diagnostics", {
  data <- make_irs_test_data()
  reversed <- data$metadata[rev(rownames(data$metadata)), , drop = FALSE]
  fit <- fit_irs(
    data$counts, reversed, "group", verbose = FALSE, seed = 1
  )
  info <- diagnostics(fit)

  expect_identical(rownames(info$samples), colnames(data$counts))
  expect_equal(info$samples$reference_total, unname(reference_totals(fit)))
  expect_length(info$change_history, fit$iterations)
})

test_that("legacy selector remains compatible with fit_irs", {
  data <- make_irs_test_data()
  selected <- select_reference_irs(
    data$counts, data$metadata, "group", verbose = FALSE, seed = 1
  )
  fit <- make_irs_test_fit(data)
  expect_identical(unname(selected), unname(fit$reference))
})

test_that("count and name validation is strict", {
  data <- make_irs_test_data()
  bad <- data$counts
  rownames(bad) <- NULL
  expect_error(
    fit_irs(bad, data$metadata, "group", verbose = FALSE),
    "taxon row names"
  )

  bad <- data$counts
  bad[1, 1] <- 1.5
  expect_error(
    fit_irs(bad, data$metadata, "group", verbose = FALSE),
    "integer counts"
  )
})
