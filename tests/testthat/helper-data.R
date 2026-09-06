make_irs_test_data <- function(n_taxa = 24L, n_samples = 50L, seed = 42L) {
  set.seed(seed)
  counts <- matrix(
    stats::rpois(n_taxa * n_samples, lambda = 25),
    nrow = n_taxa,
    ncol = n_samples,
    dimnames = list(
      paste0("taxon", seq_len(n_taxa)),
      paste0("sample", seq_len(n_samples))
    )
  )
  metadata <- data.frame(
    group = factor(
      rep(c("control", "case"), length.out = n_samples),
      levels = c("control", "case")
    ),
    row.names = colnames(counts)
  )
  list(counts = counts, metadata = metadata)
}

make_irs_test_fit <- function(data = make_irs_test_data()) {
  fit_irs(
    data$counts,
    data$metadata,
    predictor = "group",
    verbose = FALSE,
    seed = 1
  )
}
