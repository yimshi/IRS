.is_default_rownames <- function(x) {
  identical(rownames(x), as.character(seq_len(nrow(x))))
}

.validate_irs_counts <- function(counts, require_names = TRUE,
                                 integerish = TRUE, object_name = "counts") {
  if (is.data.frame(counts)) counts <- as.matrix(counts)
  if (is.null(dim(counts)) || length(dim(counts)) != 2L) {
    stop(object_name, " must be a two-dimensional count matrix.", call. = FALSE)
  }
  counts <- as.matrix(counts)
  storage.mode(counts) <- "double"
  if (!nrow(counts) || !ncol(counts)) {
    stop(object_name, " must contain at least one taxon and one sample.",
         call. = FALSE)
  }
  if (any(!is.finite(counts)) || any(counts < 0)) {
    stop(object_name, " must contain finite, non-negative values.",
         call. = FALSE)
  }
  if (integerish && any(abs(counts - round(counts)) > 1e-8)) {
    stop(object_name, " must contain integer counts.", call. = FALSE)
  }
  if (require_names) {
    if (is.null(rownames(counts)) || anyNA(rownames(counts)) ||
        any(!nzchar(rownames(counts))) || anyDuplicated(rownames(counts))) {
      stop(object_name, " must have unique, non-empty taxon row names.",
           call. = FALSE)
    }
    if (is.null(colnames(counts)) || anyNA(colnames(counts)) ||
        any(!nzchar(colnames(counts))) || anyDuplicated(colnames(counts))) {
      stop(object_name, " must have unique, non-empty sample column names.",
           call. = FALSE)
    }
  }
  counts
}

.align_irs_metadata <- function(counts, meta_dat) {
  if (!is.data.frame(meta_dat)) meta_dat <- as.data.frame(meta_dat)
  if (nrow(meta_dat) != ncol(counts)) {
    stop("meta_dat must have one row for every sample in counts.", call. = FALSE)
  }
  samples <- colnames(counts)
  if (is.null(samples)) return(meta_dat)

  if (!is.null(rownames(meta_dat)) && all(samples %in% rownames(meta_dat))) {
    return(meta_dat[samples, , drop = FALSE])
  }
  if (.is_default_rownames(meta_dat)) {
    rownames(meta_dat) <- samples
    return(meta_dat)
  }
  stop(
    "Metadata row names do not match the count-matrix sample names. ",
    "Provide matching names or metadata in the same order with default row names.",
    call. = FALSE
  )
}

.check_irs_fit <- function(fit) {
  if (!inherits(fit, "irs_fit")) {
    stop("fit must be an object returned by fit_irs().", call. = FALSE)
  }
  invisible(fit)
}

.geometric_mean <- function(x) {
  if (any(!is.finite(x)) || any(x <= 0)) {
    stop("A geometric mean requires finite positive values.", call. = FALSE)
  }
  exp(mean(log(x)))
}

.resolve_irs_counts <- function(fit, counts = NULL, integerish = TRUE) {
  .check_irs_fit(fit)
  if (is.null(counts)) counts <- fit$counts
  if (is.null(counts)) {
    stop("Supply counts, or create fit with store_counts = TRUE.", call. = FALSE)
  }
  counts <- .validate_irs_counts(
    counts, require_names = TRUE, integerish = integerish
  )
  missing_samples <- setdiff(fit$sample_names, colnames(counts))
  extra_samples <- setdiff(colnames(counts), fit$sample_names)
  if (length(missing_samples) || length(extra_samples)) {
    stop(
      "Count-matrix samples do not match the IRS fit. Missing: ",
      paste(missing_samples, collapse = ", "), "; extra: ",
      paste(extra_samples, collapse = ", "), ".",
      call. = FALSE
    )
  }
  counts[, fit$sample_names, drop = FALSE]
}

#' Fit IRS normalization
#'
#' Runs iterative reference selection and stores the final reference taxa,
#' reference totals, centered size factors, offsets, convergence information,
#' and sample diagnostics needed by downstream adapters.
#'
#' @param counts Numeric count matrix with taxa in rows and samples in columns.
#'   Unique taxon and sample names are required.
#' @param meta_dat Sample metadata. Matching row names are preferred; metadata
#'   with default row names is assumed to already follow the count-column order.
#' @param predictor Name of the primary predictor in `meta_dat`.
#' @param max_iter,tolerance,quantile_range,p_threshold,convergence_window
#'   Arguments passed to [select_reference_irs()].
#' @param verbose,seed,kendall_transform,initial_ref Arguments passed to
#'   [select_reference_irs()].
#' @param store_counts Store a copy of the count matrix in the fitted object.
#'   The default is `FALSE` to keep large fits small.
#'
#' @return An object of class `irs_fit`.
#' @export
fit_irs <- function(
    counts,
    meta_dat,
    predictor,
    max_iter = 30,
    tolerance = 1e-3,
    quantile_range = c(0, 1),
    p_threshold = 0.05,
    convergence_window = 3,
    verbose = TRUE,
    seed = NULL,
    kendall_transform = c("none", "log1p"),
    initial_ref = NULL,
    store_counts = FALSE) {
  counts <- .validate_irs_counts(counts, require_names = TRUE, integerish = TRUE)
  if (nrow(counts) < 2L || ncol(counts) < 2L) {
    stop("IRS requires at least two taxa and two samples.", call. = FALSE)
  }
  meta_dat <- .align_irs_metadata(counts, meta_dat)
  kendall_transform <- match.arg(kendall_transform)

  details <- .select_reference_irs_details(
    otu_tab_sim = counts,
    meta_dat = meta_dat,
    predictor = predictor,
    max_iter = max_iter,
    tolerance = tolerance,
    quantile_range = quantile_range,
    p_threshold = p_threshold,
    convergence_window = convergence_window,
    verbose = verbose,
    seed = seed,
    kendall_transform = kendall_transform,
    initial_ref = initial_ref
  )
  reference <- details$reference
  names(reference) <- rownames(counts)
  reference_total <- colSums(counts[reference, , drop = FALSE])
  if (any(reference_total <= 0)) {
    bad <- names(reference_total)[reference_total <= 0]
    stop(
      "The final IRS reference total is zero for sample(s): ",
      paste(bad, collapse = ", "),
      ". IRS normalization requires a positive reference total.",
      call. = FALSE
    )
  }
  size_factor <- reference_total / .geometric_mean(reference_total)
  library_total <- colSums(counts)

  fit <- list(
    call = match.call(),
    predictor = predictor,
    reference = reference,
    reference_taxa = names(reference)[reference],
    reference_total = reference_total,
    size_factors = size_factor,
    log_offsets = log(size_factor),
    library_total = library_total,
    reference_fraction = ifelse(
      library_total > 0, reference_total / library_total, NA_real_
    ),
    converged = details$converged,
    iterations = details$iterations,
    change_history = details$change_history,
    reference_size_history = details$reference_size_history,
    feature_names = rownames(counts),
    sample_names = colnames(counts),
    count_dim = dim(counts),
    counts = if (isTRUE(store_counts)) counts else NULL,
    package_version = utils::packageVersion("IRS")
  )
  class(fit) <- "irs_fit"
  if (!fit$converged) {
    warning(
      "IRS reached max_iter without satisfying the convergence criterion.",
      call. = FALSE
    )
  }
  fit
}

#' @export
print.irs_fit <- function(x, ...) {
  cat("IRS fit\n")
  cat("  taxa:", x$count_dim[1], "\n")
  cat("  samples:", x$count_dim[2], "\n")
  cat("  reference taxa:", length(x$reference_taxa), "\n")
  cat("  iterations:", x$iterations, "\n")
  cat("  converged:", x$converged, "\n")
  invisible(x)
}

#' IRS fit accessors
#'
#' @param fit An `irs_fit` object.
#' @param counts Optional count matrix. Required when counts were not stored in
#'   `fit`.
#' @param orientation Orientation of the returned normalized count matrix.
#' @param taxa_are_rows Whether taxa occur in rows of `counts`.
#'
#' @return `reference_taxa()` returns character names; `reference_indices()`
#'   returns integer positions; `reference_totals()`, `size_factors()`, and
#'   `offsets()` return named numeric vectors; `normalized_counts()` returns a
#'   numeric matrix; `diagnostics()` returns a diagnostic list.
#' @name irs_accessors
NULL

#' @rdname irs_accessors
#' @export
reference_taxa <- function(fit) {
  .check_irs_fit(fit)
  fit$reference_taxa
}

#' @rdname irs_accessors
#' @export
reference_indices <- function(fit, counts, taxa_are_rows = TRUE) {
  .check_irs_fit(fit)
  counts <- as.matrix(counts)
  taxa_names <- if (isTRUE(taxa_are_rows)) rownames(counts) else colnames(counts)
  if (is.null(taxa_names) || anyDuplicated(taxa_names)) {
    stop("The taxon dimension must have unique names.", call. = FALSE)
  }
  missing <- setdiff(fit$reference_taxa, taxa_names)
  if (length(missing)) {
    stop("Reference taxa missing from counts: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  match(fit$reference_taxa, taxa_names)
}

#' @rdname irs_accessors
#' @export
reference_totals <- function(fit) {
  .check_irs_fit(fit)
  fit$reference_total
}

#' @rdname irs_accessors
#' @export
size_factors <- function(fit) {
  .check_irs_fit(fit)
  fit$size_factors
}

#' @rdname irs_accessors
#' @export
offsets <- function(fit) {
  .check_irs_fit(fit)
  fit$log_offsets
}

#' @rdname irs_accessors
#' @export
normalized_counts <- function(
    fit, counts = NULL,
    orientation = c("taxa_by_samples", "samples_by_taxa")) {
  orientation <- match.arg(orientation)
  counts <- .resolve_irs_counts(fit, counts, integerish = TRUE)
  out <- sweep(counts, 2L, fit$size_factors, "/")
  if (orientation == "samples_by_taxa") out <- t(out)
  out
}

#' @rdname irs_accessors
#' @export
diagnostics <- function(fit) {
  .check_irs_fit(fit)
  list(
    converged = fit$converged,
    iterations = fit$iterations,
    change_history = fit$change_history,
    reference_size_history = fit$reference_size_history,
    samples = data.frame(
      sample = fit$sample_names,
      library_total = unname(fit$library_total),
      reference_total = unname(fit$reference_total),
      reference_fraction = unname(fit$reference_fraction),
      size_factor = unname(fit$size_factors),
      row.names = fit$sample_names,
      check.names = FALSE
    )
  )
}
