#' Iterative Reference Selection (IRS) for Microbiome Data
#'
#' @description
#' This function implements the Iterative Reference Selection (IRS) algorithm.
#' It identifies a stable set of reference taxa that are not associated with
#' the predictor of interest, which can then be used for normalization in
#' downstream differential abundance analysis.
#'
#' @details
#' The selection process is an iterative procedure that combines two robust
#' statistical filters:
#' \itemize{
#'   \item \strong{Kendall's Rank Correlation}: A non-parametric test to assess
#'   monotonic relationships on normalized abundances.
#'   \item \strong{Robust Poisson GLM}: A Generalized Linear Model with sandwich
#'   (HC3) variance estimators to account for overdispersion and potential
#'   model misspecification in count data.
#' }
#' Convergence is achieved when the proportion of change in the reference set
#' falls below the specified \code{tolerance}.
#'
#' @param otu_tab_sim A numeric \code{matrix} of taxon counts (taxa as rows,
#' samples as columns).
#' @param meta_dat A \code{data.frame} containing sample metadata.
#' @param predictor A \code{string} specifying the column name in \code{meta_dat}
#' to be used as the primary predictor (exposure).
#' @param max_iter Integer. Maximum number of iterations allowed. Default is 30.
#' @param tolerance Numeric. Convergence threshold for the average change in the
#' reference set. Default is 1e-3.
#' @param quantile_range A numeric vector of length 2. Defines the abundance
#' quantile range for the initial reference set. Default is \code{c(0, 1)}.
#' @param p_threshold Numeric. The p-value threshold for both Kendall and GLM
#' filters. Taxa with p-values above this threshold are considered "stable".
#' Default is 0.05.
#' @param convergence_window Integer. Number of recent iterations to average for
#' checking convergence. Default is 3.
#' @param verbose Logical. If \code{TRUE}, prints iteration progress to the
#' console.
#' @param seed Integer. Random seed for reproducibility.
#' @param kendall_transform Character. Transformation to apply before Kendall
#' test: \code{"none"} or \code{"log1p"}.
#' @param initial_ref Optional. A logical vector or character vector of taxa
#' names to seed the first iteration.
#'
#' @return A logical vector of length \code{nrow(otu_tab_sim)} indicating
#' whether each taxon belongs to the final reference set.
#'
#' @importFrom stats model.matrix reformulate quantile cor.test pnorm vcov glm poisson var
#' @importFrom fastglm fastglm
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming otu_matrix and metadata are pre-loaded
#' ref_set <- select_reference_irs(
#'   otu_tab_sim = otu_matrix,
#'   meta_dat = metadata,
#'   predictor = "Group",
#'   p_threshold = 0.05
#' )
#' }

select_reference_irs <- function(
    otu_tab_sim,           # matrix: taxa x samples (counts)
    meta_dat,              # data.frame: sample metadata
    predictor,             # string: column name in meta_dat
    max_iter = 30,
    tolerance = 1e-3,
    quantile_range = c(0, 1),
    p_threshold = 0.05,    # final Kendall∧GLM threshold
    convergence_window = 3,
    verbose = TRUE,
    seed = NULL,
    kendall_transform = c("none","log1p"),
    initial_ref = NULL     # NULL | logical(length = n_taxa) | character(rownames)
) {
  if (!is.null(seed)) set.seed(seed)
  kendall_transform <- match.arg(kendall_transform)

  # basic shapes
  if (is.null(dim(otu_tab_sim))) stop("otu_tab_sim must be a matrix (taxa x samples).")
  n_taxa  <- nrow(otu_tab_sim)
  n_samps <- ncol(otu_tab_sim)
  if (n_taxa < 2L || n_samps < 2L) stop("Need at least 2 taxa and 2 samples.")

  # design matrix
  if (!is.data.frame(meta_dat)) meta_dat <- as.data.frame(meta_dat)
  if (!(predictor %in% colnames(meta_dat))) {
    stop(sprintf("predictor '%s' not found in meta_dat.", predictor))
  }
  X <- model.matrix(reformulate(predictor, response = NULL), data = meta_dat)

  # seed reference set
  taxa_abundance <- rowSums(otu_tab_sim)
  if (is.null(initial_ref)) {
    q_bounds <- stats::quantile(taxa_abundance, probs = quantile_range)
    ref_taxa <- (taxa_abundance >= q_bounds[1]) & (taxa_abundance <= q_bounds[2])
  } else if (is.logical(initial_ref)) {
    if (length(initial_ref) != n_taxa) {
      stop("initial_ref logical vector must have length equal to nrow(otu_tab_sim).")
    }
    ref_taxa <- initial_ref
  } else if (is.character(initial_ref)) {
    if (is.null(rownames(otu_tab_sim))) {
      stop("Row names are required in otu_tab_sim when initial_ref is character.")
    }
    ref_taxa <- rownames(otu_tab_sim) %in% initial_ref
  } else {
    stop("initial_ref must be NULL, logical, or character.")
  }
  if (!any(ref_taxa)) stop("Initial reference set is empty after seeding.")

  change_history <- numeric(0)
  iter <- 0L

  # one iteration: recompute offsets, p-values, and update reference set
  step_once <- function(threshold_p) {
    total_ref <- colSums(otu_tab_sim[ref_taxa, , drop = FALSE])
    denom_log <- precompute_log_offsets(total_ref, otu_tab_sim, ref_taxa, eps = 1e-8)

    # Kendall on normalized abundances (computed per-taxon to avoid large temp matrix)
    p_rank <- vapply(
      seq_len(n_taxa),
      FUN = function(i) {
        y_norm <- otu_tab_sim[i, ] / exp(denom_log[i, ])
        p <- kendall_pval(y_norm, X, transform = kendall_transform)
        if (!is.finite(p)) 1 else p
      },
      FUN.VALUE = numeric(1)
    )

    # Poisson-GLM with sandwich variance on raw counts + offsets
    p_glm <- vapply(
      seq_len(n_taxa),
      FUN = function(i) {
        p <- robust_pval(otu_tab_sim[i, ], X, denom_log[i, ])
        if (!is.finite(p)) 1 else p
      },
      FUN.VALUE = numeric(1)
    )

    new_ref <- (p_rank >= threshold_p) & (p_glm >= threshold_p)
    list(new_ref = new_ref)
  }

  # main convergence loop (single threshold = p_threshold)
  for (k in seq_len(max_iter)) {
    iter <- iter + 1L
    st <- step_once(threshold_p = p_threshold)
    new_ref_taxa <- st$new_ref

    delta <- mean(new_ref_taxa != ref_taxa)
    change_history <- c(change_history, delta)
    if (verbose) cat("Iter:", k, "Reference size:", sum(new_ref_taxa),
                     "Δ:", round(delta, 4), "\n")

    # update, then check averaged recent change
    ref_taxa <- new_ref_taxa
    if (length(change_history) >= convergence_window) {
      avg_delta <- mean(tail(change_history, convergence_window))
      if (avg_delta < tolerance) {
        if (verbose) message(sprintf("Converged at iteration %d", iter))
        break
      }
    }
  }

  return(ref_taxa)
}

#' @keywords internal
kendall_pval <- function(
    y,                 # numeric: normalized abundance for one taxon (length = n samples)
    X,                 # design matrix; use column 2 as predictor by default
    offset_vec = NULL, # unused; kept for API compatibility
    x_col = 2,
    transform = c("none","log1p")  # optional transform to reduce ties/heavy tails
) {
  # --- helpers ---------------------------------------------------------------
  # Kendall p-value helper (robust to extremes; ties handled internally)
  transform <- match.arg(transform)
  # predictor vector
  stopifnot(is.matrix(X) || is.data.frame(X))
  if (x_col > ncol(X)) stop("x_col exceeds number of columns in X")
  x <- X[, x_col]
  if (is.factor(x))    x <- as.numeric(as.character(x))
  if (is.character(x)) x <- as.numeric(x)

  # optional transform on response
  if (transform == "log1p") y <- log1p(y)

  # finite pairs only
  keep <- is.finite(y) & is.finite(x)
  y <- y[keep]; x <- x[keep]

  # need variability
  if (length(y) < 3L || length(unique(y)) < 2L || length(unique(x)) < 2L) {
    return(NA_real_)
  }
  suppressWarnings(as.numeric(stats::cor.test(x, y, method = "kendall", exact = FALSE)$p.value))
}

#' @keywords internal
robust_pval <- function(y, X, offset_vec,
                        mu_floor = 1e-10,
                        ridge    = 1e-8) {
  # Poisson-GLM with sandwich (HC3) p-value for the 2nd coefficient
  if (length(y) < 3L || var(y) == 0 || ncol(X) < 2L) return(1)
  fit <- tryCatch(
    fastglm::fastglm(x = X, y = y, family = poisson(),
                     offset = offset_vec, model = FALSE),
    error = function(e) NULL
  )
  if (is.null(fit)) return(1)

  beta_hat <- fit$coefficients
  if (!all(is.finite(beta_hat))) return(1)

  eta <- offset_vec + as.vector(X %*% beta_hat)
  mu  <- exp(eta)
  mu  <- pmax(mu, mu_floor)
  WX   <- X * mu
  XtWX <- crossprod(X, WX)
  p <- ncol(X)
  XtWX_reg <- XtWX + diag(ridge, p)

  bread <- tryCatch(solve(XtWX_reg),
                    error = function(e) NULL)
  if (is.null(bread) || any(!is.finite(bread))) {
    V_fisher <- tryCatch(stats::vcov(glm(y ~ X[,2], family = poisson(),
                                         offset = offset_vec)),
                         error = function(e) NULL)
    if (is.null(V_fisher)) return(0)
    se <- sqrt(diag(V_fisher))
    if (!is.finite(se[2]) || se[2] == 0) return(0)
    z <- beta_hat[2] / se[2]
    return(2 * stats::pnorm(-abs(z)))
  }

  # Meat: X' diag((y-mu)^2) X
  XR   <- X * (y - mu)
  meat <- crossprod(XR)

  vcov <- bread %*% meat %*% bread
  if (any(!is.finite(vcov))) return(0)

  se2 <- diag(vcov)
  if (length(se2) < 2L || !is.finite(se2[2]) || se2[2] <= 0) return(0)

  z <- beta_hat[2] / sqrt(se2[2])
  2 * stats::pnorm(-abs(z))
}

#' @keywords internal
precompute_log_offsets <- function(total_ref, otu_tab, ref_taxa, eps = 1e-8) {
  # precompute log-offset matrix: log(total_ref - self), with floor eps
  off <- matrix(total_ref, nrow = nrow(otu_tab), ncol = ncol(otu_tab), byrow = TRUE)
  if (any(ref_taxa)) {
    off[ref_taxa, ] <- off[ref_taxa, ] - otu_tab[ref_taxa, , drop = FALSE]
  }
  log(pmax(off, eps))
}

