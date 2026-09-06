#' Iterative Reference Selection (IRS) for Microbiome Data
#'
#' Identifies a stable set of reference taxa that are not associated with the
#' predictor of interest. The returned logical vector can be used directly by
#' existing code; [fit_irs()] provides the richer interface used by downstream
#' adapters.
#'
#' @param otu_tab_sim Numeric matrix of counts with taxa in rows and samples in
#'   columns.
#' @param meta_dat Data frame containing one row per sample.
#' @param predictor Name of the primary predictor in `meta_dat`.
#' @param max_iter Maximum number of iterations.
#' @param tolerance Convergence threshold for the recent mean proportion of
#'   reference-set changes.
#' @param quantile_range Two probabilities defining the abundance range used to
#'   seed the initial reference set.
#' @param p_threshold Taxa with both IRS screening p-values greater than or
#'   equal to this value remain in the reference set.
#' @param convergence_window Number of recent iterations used for convergence.
#' @param verbose Print iteration progress.
#' @param seed Optional random seed retained for reproducible workflows.
#' @param kendall_transform Transformation used before the Kendall screen.
#' @param initial_ref Optional logical vector or character vector of taxa names
#'   used to seed the first iteration.
#'
#' @return A logical vector indicating membership in the final reference set.
#' @export
#' @importFrom fastglm fastglm
#' @importFrom stats cor.test glm model.matrix pnorm poisson quantile
#'   reformulate var vcov
#'
#' @examples
#' \dontrun{
#' ref <- select_reference_irs(counts, metadata, predictor = "group")
#' }
select_reference_irs <- function(
    otu_tab_sim,
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
    initial_ref = NULL) {
  details <- .select_reference_irs_details(
    otu_tab_sim = otu_tab_sim,
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
  details$reference
}

.select_reference_irs_details <- function(
    otu_tab_sim,
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
    initial_ref = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  kendall_transform <- match.arg(kendall_transform)

  if (is.null(dim(otu_tab_sim))) {
    stop("otu_tab_sim must be a matrix with taxa in rows and samples in columns.",
         call. = FALSE)
  }
  otu_tab_sim <- as.matrix(otu_tab_sim)
  storage.mode(otu_tab_sim) <- "double"
  n_taxa <- nrow(otu_tab_sim)
  n_samps <- ncol(otu_tab_sim)
  if (n_taxa < 2L || n_samps < 2L) {
    stop("IRS requires at least two taxa and two samples.", call. = FALSE)
  }
  if (any(!is.finite(otu_tab_sim)) || any(otu_tab_sim < 0)) {
    stop("otu_tab_sim must contain finite, non-negative counts.", call. = FALSE)
  }

  if (!is.data.frame(meta_dat)) {
    meta_dat <- as.data.frame(meta_dat)
  }
  if (nrow(meta_dat) != n_samps) {
    stop("meta_dat must have one row per sample in otu_tab_sim.", call. = FALSE)
  }
  if (length(predictor) != 1L || !predictor %in% colnames(meta_dat)) {
    stop(sprintf("predictor '%s' was not found in meta_dat.", predictor),
         call. = FALSE)
  }
  X <- stats::model.matrix(stats::reformulate(predictor), data = meta_dat)
  if (nrow(X) != n_samps) {
    stop("Missing predictor values are not supported; remove or impute them first.",
         call. = FALSE)
  }
  if (ncol(X) < 2L) {
    stop("The predictor does not generate a testable model coefficient.",
         call. = FALSE)
  }

  if (length(max_iter) != 1L || max_iter < 1 || max_iter != as.integer(max_iter)) {
    stop("max_iter must be a positive integer.", call. = FALSE)
  }
  if (length(convergence_window) != 1L || convergence_window < 1 ||
      convergence_window != as.integer(convergence_window)) {
    stop("convergence_window must be a positive integer.", call. = FALSE)
  }
  if (length(tolerance) != 1L || !is.finite(tolerance) || tolerance < 0) {
    stop("tolerance must be a finite non-negative number.", call. = FALSE)
  }
  if (length(p_threshold) != 1L || !is.finite(p_threshold) ||
      p_threshold < 0 || p_threshold > 1) {
    stop("p_threshold must lie between zero and one.", call. = FALSE)
  }
  if (length(quantile_range) != 2L || any(!is.finite(quantile_range)) ||
      any(quantile_range < 0 | quantile_range > 1) ||
      quantile_range[1] > quantile_range[2]) {
    stop("quantile_range must contain two ordered probabilities.", call. = FALSE)
  }

  taxa_abundance <- rowSums(otu_tab_sim)
  if (is.null(initial_ref)) {
    q_bounds <- stats::quantile(taxa_abundance, probs = quantile_range)
    ref_taxa <- taxa_abundance >= q_bounds[1] & taxa_abundance <= q_bounds[2]
  } else if (is.logical(initial_ref)) {
    if (length(initial_ref) != n_taxa || anyNA(initial_ref)) {
      stop("initial_ref must be a non-missing logical vector with one value per taxon.",
           call. = FALSE)
    }
    ref_taxa <- initial_ref
  } else if (is.character(initial_ref)) {
    if (is.null(rownames(otu_tab_sim))) {
      stop("Taxon row names are required when initial_ref contains names.",
           call. = FALSE)
    }
    unknown <- setdiff(initial_ref, rownames(otu_tab_sim))
    if (length(unknown)) {
      stop("Unknown taxa in initial_ref: ", paste(unknown, collapse = ", "),
           call. = FALSE)
    }
    ref_taxa <- rownames(otu_tab_sim) %in% initial_ref
  } else {
    stop("initial_ref must be NULL, logical, or character.", call. = FALSE)
  }
  if (!any(ref_taxa)) {
    stop("The initial reference set is empty.", call. = FALSE)
  }

  change_history <- numeric()
  reference_size_history <- integer()
  converged <- FALSE

  for (iteration in seq_len(max_iter)) {
    total_ref <- colSums(otu_tab_sim[ref_taxa, , drop = FALSE])
    denom_log <- precompute_log_offsets(
      total_ref, otu_tab_sim, ref_taxa, eps = 1e-8
    )

    p_rank <- vapply(seq_len(n_taxa), function(i) {
      y_norm <- otu_tab_sim[i, ] / exp(denom_log[i, ])
      p <- kendall_pval(y_norm, X, transform = kendall_transform)
      if (is.finite(p)) p else 1
    }, numeric(1))

    p_glm <- vapply(seq_len(n_taxa), function(i) {
      p <- robust_pval(otu_tab_sim[i, ], X, denom_log[i, ])
      if (is.finite(p)) p else 1
    }, numeric(1))

    new_ref_taxa <- p_rank >= p_threshold & p_glm >= p_threshold
    if (!any(new_ref_taxa)) {
      stop("IRS produced an empty reference set at iteration ", iteration, ".",
           call. = FALSE)
    }

    delta <- mean(new_ref_taxa != ref_taxa)
    change_history <- c(change_history, delta)
    reference_size_history <- c(reference_size_history, sum(new_ref_taxa))
    if (verbose) {
      cat("Iter:", iteration, "Reference size:", sum(new_ref_taxa),
          "change:", round(delta, 4), "\n")
    }

    ref_taxa <- new_ref_taxa
    if (length(change_history) >= convergence_window &&
        mean(utils::tail(change_history, convergence_window)) < tolerance) {
      converged <- TRUE
      if (verbose) {
        message(sprintf("Converged at iteration %d", iteration))
      }
      break
    }
  }

  names(ref_taxa) <- rownames(otu_tab_sim)
  list(
    reference = ref_taxa,
    converged = converged,
    iterations = length(change_history),
    change_history = change_history,
    reference_size_history = reference_size_history
  )
}

#' @keywords internal
kendall_pval <- function(y, X, offset_vec = NULL, x_col = 2,
                         transform = c("none", "log1p")) {
  transform <- match.arg(transform)
  stopifnot(is.matrix(X) || is.data.frame(X))
  if (x_col > ncol(X)) {
    stop("x_col exceeds the number of columns in X.", call. = FALSE)
  }
  x <- X[, x_col]
  if (is.factor(x)) x <- as.numeric(as.character(x))
  if (is.character(x)) x <- as.numeric(x)
  if (transform == "log1p") y <- log1p(y)

  keep <- is.finite(y) & is.finite(x)
  y <- y[keep]
  x <- x[keep]
  if (length(y) < 3L || length(unique(y)) < 2L || length(unique(x)) < 2L) {
    return(NA_real_)
  }
  suppressWarnings(as.numeric(stats::cor.test(
    x, y, method = "kendall", exact = FALSE
  )$p.value))
}

#' @keywords internal
robust_pval <- function(y, X, offset_vec, mu_floor = 1e-10, ridge = 1e-8) {
  if (length(y) < 3L || stats::var(y) == 0 || ncol(X) < 2L) return(1)
  fit <- tryCatch(
    fastglm::fastglm(
      x = X, y = y, family = stats::poisson(), offset = offset_vec,
      model = FALSE
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) return(1)

  beta_hat <- fit$coefficients
  if (!all(is.finite(beta_hat))) return(1)

  eta <- offset_vec + as.vector(X %*% beta_hat)
  mu <- pmax(exp(eta), mu_floor)
  XtWX <- crossprod(X, X * mu)
  p <- ncol(X)
  bread <- tryCatch(solve(XtWX + diag(ridge, p)), error = function(e) NULL)
  if (is.null(bread) || any(!is.finite(bread))) {
    V_fisher <- tryCatch(
      stats::vcov(stats::glm(
        y ~ X[, 2], family = stats::poisson(), offset = offset_vec
      )),
      error = function(e) NULL
    )
    if (is.null(V_fisher)) return(0)
    se <- sqrt(diag(V_fisher))
    if (!is.finite(se[2]) || se[2] == 0) return(0)
    return(2 * stats::pnorm(-abs(beta_hat[2] / se[2])))
  }

  meat <- crossprod(X * (y - mu))
  vcov <- bread %*% meat %*% bread
  if (any(!is.finite(vcov))) return(0)
  se2 <- diag(vcov)
  if (length(se2) < 2L || !is.finite(se2[2]) || se2[2] <= 0) return(0)
  2 * stats::pnorm(-abs(beta_hat[2] / sqrt(se2[2])))
}

#' @keywords internal
precompute_log_offsets <- function(total_ref, otu_tab, ref_taxa, eps = 1e-8) {
  off <- matrix(
    total_ref, nrow = nrow(otu_tab), ncol = ncol(otu_tab), byrow = TRUE
  )
  if (any(ref_taxa)) {
    off[ref_taxa, ] <- off[ref_taxa, ] - otu_tab[ref_taxa, , drop = FALSE]
  }
  log(pmax(off, eps))
}
