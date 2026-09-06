# The implementation in this file is derived from ZicoSeq in GUniFrac 1.9,
# Copyright (C) Jun Chen, Xianyang Zhang, Lu Yang, and Lujun Zhang.
# GUniFrac is distributed under GPL-3. IRS modifications replace ZicoSeq's
# internal reference selection with a fixed IRS reference set, align references
# by taxon name, retain raw IRS reference totals independently of feature
# filtering/winsorization, and add explicit zero-denominator handling.

.irs_perm_fdr_adj <- function(F0, Fp) {
  ord <- order(F0, decreasing = TRUE)
  F0 <- F0[ord]
  perm_no <- ncol(Fp)
  Fp <- sort(c(as.vector(Fp)[!is.na(as.vector(Fp))], F0), decreasing = FALSE)
  n <- length(Fp)
  m <- length(F0)
  false_positive_number <- (n + 1) - match(F0, Fp) - seq_len(m)
  adjusted <- false_positive_number / perm_no / seq_len(m)
  zero <- adjusted == 0
  adjusted[zero] <- 0.5 / perm_no / which(zero)
  pmin(1, rev(cummin(rev(adjusted))))[order(ord)]
}

.irs_perm_fwer_adj <- function(F0, Fp) {
  ord <- order(F0, decreasing = TRUE)
  m <- length(F0)
  F0 <- F0[ord]
  Fp <- Fp[ord, , drop = FALSE]
  column_max <- Fp[m, ]
  adjusted <- vapply(m:1, function(i) {
    column_max <<- pmax(Fp[i, ], column_max)
    mean(c(1, column_max >= F0[i]))
  }, numeric(1))
  adjusted <- rev(adjusted)
  pmin(1, rev(cummin(rev(adjusted))))[order(ord)]
}

.irs_na_pad <- function(x, missing) {
  out <- numeric(length(missing))
  out[!missing] <- x
  out[missing] <- NA_real_
  out
}

.irs_bbmix_fit_mm <- function(count, depth, n_iter = 10, winsor_quantile = 1) {
  if (mean(count == 0) >= 0.95 || sum(count != 0) < 5) {
    stop(
      "Too few nonzero observations to fit the beta-binomial mixture.",
      call. = FALSE
    )
  }
  proportion <- count / depth
  upper <- stats::quantile(proportion, winsor_quantile)
  proportion[proportion >= upper] <- upper
  var1 <- stats::var(proportion) / 2
  mean1 <- mean(proportion) / 2
  var2 <- stats::var(proportion) / 2
  mean2 <- 3 * mean(proportion) / 2
  mixture_probability <- 0.5

  for (iteration in seq_len(n_iter)) {
    shape1.1 <- ((1 - mean1) / var1 - 1 / mean1) * mean1^2
    shape1.2 <- shape1.1 * (1 / mean1 - 1)
    shape2.1 <- ((1 - mean2) / var2 - 1 / mean2) * mean2^2
    shape2.2 <- shape2.1 * (1 / mean2 - 1)
    m1 <- shape1.1 / (shape1.1 + shape1.2)
    s1 <- shape1.1 + shape1.2
    m2 <- shape2.1 / (shape2.1 + shape2.2)
    s2 <- shape2.1 + shape2.2
    f1 <- rmutil::dbetabinom(count, depth, m1, s1)
    f2 <- rmutil::dbetabinom(count, depth, m2, s2)
    q1 <- mixture_probability * f1 /
      (mixture_probability * f1 + (1 - mixture_probability) * f2)
    q2 <- 1 - q1
    mixture_probability <- mean(q1)
    mean1 <- sum(proportion * q1) / sum(q1)
    var1 <- sum((proportion - mean1)^2 * q1) / sum(q1)
    mean2 <- sum(proportion * q2) / sum(q2)
    var2 <- sum((proportion - mean2)^2 * q2) / sum(q2)
  }
  list(
    shape1.1 = shape1.1,
    shape1.2 = shape1.2,
    shape2.1 = shape2.1,
    shape2.2 = shape2.2,
    pi = mixture_probability,
    q1 = q1
  )
}

#' ZicoSeq with a fixed IRS reference set
#'
#' Runs the ZicoSeq count-data pipeline while replacing its internal reference
#' selection with taxa selected by IRS. Reference taxa are aligned by name and
#' their totals are computed from the original input counts before ZicoSeq
#' feature filtering or winsorization. The reference set is fixed for the
#' single testing stage; posterior sampling, covariate adjustment, link
#' functions, stratified permutation, and ZicoSeq-style multiplicity
#' adjustment are retained.
#'
#' Supply the original integer count matrix rather than output from
#' [normalized_counts()]. Use this wrapper in place of a subsequent call to
#' `GUniFrac::ZicoSeq()`: native ZicoSeq reference-selection controls are not
#' used because IRS supplies a fixed set. The original-count reference totals
#' remain fixed even when ZicoSeq later filters or winsorizes tested features.
#' The default `zero_ref_action = "error"` preserves the IRS denominator
#' definition; the fallback option changes zero denominators and is intended
#' only for an explicitly reported sensitivity analysis.
#'
#' This function is derived from `GUniFrac::ZicoSeq()` version 1.9 under GPL-3.
#'
#' @param meta.dat Sample metadata with samples in rows.
#' @param feature.dat Integer count matrix with taxa in rows and samples in
#'   columns.
#' @param grp.name Name of the variable of interest in `meta.dat`.
#' @param adj.name Optional adjustment-variable names.
#' @param feature.dat.type Retained for interface familiarity. IRS integration
#'   currently supports count data only.
#' @param prev.filter,mean.abund.filter,max.abund.filter ZicoSeq feature filters.
#' @param min.prop Lower bound for sampled proportions.
#' @param is.winsor,outlier.pct,winsor.end ZicoSeq winsorization settings.
#' @param is.post.sample,post.sample.no Posterior-sampling settings.
#' @param link.func Link functions applied before testing.
#' @param stats.combine.func Function used to combine link-function statistics.
#' @param perm.no Number of permutations.
#' @param strata Optional permutation strata or a metadata column name.
#' @param ref.pct Retained for compatibility and ignored because IRS supplies
#'   fixed references.
#' @param stage.no Retained for compatibility. Fixed IRS references require one
#'   stage and values other than one are ignored with a warning.
#' @param excl.pct Retained for compatibility and ignored.
#' @param p.max Retained for compatibility and ignored.
#' @param is.fwer Compute permutation FWER-adjusted p-values.
#' @param verbose Print progress.
#' @param return.feature.dat Return the processed feature matrix.
#' @param irs_fit An optional `irs_fit` object. Supply exactly one of `irs_fit`
#'   and `ref.features`.
#' @param ref.features Optional character vector of IRS reference-taxon names.
#' @param zero_ref_action Action when an original IRS reference total is zero.
#'   The default stops; `"zicoseq_fallback"` uses half the smallest positive
#'   divisor and records the affected samples.
#'
#' @return A list compatible with the principal fields returned by ZicoSeq,
#'   plus an `irs` provenance and normalization component.
#' @export
#' @import vegan
#' @importFrom matrixStats rowMaxs rowSds
zicoseq_irs <- function(
    meta.dat,
    feature.dat,
    grp.name,
    adj.name = NULL,
    feature.dat.type = c("count", "proportion", "other"),
    prev.filter = 0,
    mean.abund.filter = 0,
    max.abund.filter = 0,
    min.prop = 0,
    is.winsor = FALSE,
    outlier.pct = 0.03,
    winsor.end = c("top", "bottom", "both"),
    is.post.sample = TRUE,
    post.sample.no = 25,
    link.func = list(function(x) sign(x) * abs(x)^0.5),
    stats.combine.func = max,
    perm.no = 99,
    strata = NULL,
    ref.pct = 0.5,
    stage.no = 1,
    excl.pct = 0,
    p.max = 500,
    is.fwer = FALSE,
    verbose = TRUE,
    return.feature.dat = TRUE,
    irs_fit = NULL,
    ref.features = NULL,
    zero_ref_action = c("error", "zicoseq_fallback")) {
  this_call <- match.call()
  feature.dat.type <- match.arg(feature.dat.type)
  winsor.end <- match.arg(winsor.end)
  zero_ref_action <- match.arg(zero_ref_action)
  if (feature.dat.type != "count") {
    stop("zicoseq_irs() currently supports feature.dat.type = 'count' only.",
         call. = FALSE)
  }
  if (!identical(stage.no, 1) || !identical(excl.pct, 0)) {
    warning(
      "Fixed IRS references use stage.no = 1 and excl.pct = 0; supplied ",
      "values were ignored.", call. = FALSE
    )
  }
  stage.no <- 1L
  excl.pct <- 0
  invisible(ref.pct)
  invisible(p.max)
  invisible(stage.no)
  invisible(excl.pct)

  feature.dat <- .validate_irs_counts(
    feature.dat, require_names = TRUE, integerish = TRUE,
    object_name = "feature.dat"
  )
  meta.dat <- .align_irs_metadata(feature.dat, meta.dat)
  if (length(grp.name) != 1L || !grp.name %in% colnames(meta.dat)) {
    stop("grp.name was not found in meta.dat.", call. = FALSE)
  }
  if (!is.null(adj.name) && !all(adj.name %in% colnames(meta.dat))) {
    stop("One or more adj.name variables were not found in meta.dat.",
         call. = FALSE)
  }
  if (xor(is.null(irs_fit), is.null(ref.features)) == FALSE) {
    stop("Supply exactly one of irs_fit and ref.features.", call. = FALSE)
  }
  if (!is.null(irs_fit)) {
    .check_irs_fit(irs_fit)
    .aligned_fit_vector(
      irs_fit, colnames(feature.dat), irs_fit$reference_total, "feature.dat"
    )
    ref.features <- irs_fit$reference_taxa
  }
  ref.features <- unique(as.character(ref.features))
  if (!length(ref.features) || anyNA(ref.features) || any(!nzchar(ref.features))) {
    stop("ref.features must contain at least one non-empty taxon name.",
         call. = FALSE)
  }
  missing_reference <- setdiff(ref.features, rownames(feature.dat))
  if (length(missing_reference)) {
    stop(
      "IRS reference taxa missing from feature.dat: ",
      paste(missing_reference, collapse = ", "), ".", call. = FALSE
    )
  }
  raw_reference_total <- colSums(
    feature.dat[ref.features, , drop = FALSE]
  )

  if (!is.null(strata)) {
    if (length(strata) == 1L && is.character(strata)) {
      if (!strata %in% colnames(meta.dat)) {
        stop("The strata column was not found in meta.dat.", call. = FALSE)
      }
      strata <- factor(meta.dat[[strata]])
    } else {
      strata <- factor(strata)
    }
    if (length(strata) != nrow(meta.dat)) {
      stop("strata must have one value per sample.", call. = FALSE)
    }
  }

  feature_sd <- matrixStats::rowSds(feature.dat)
  if (any(feature_sd == 0)) {
    stop(
      "Features with identical values must be removed before ZicoSeq: ",
      paste(rownames(feature.dat)[feature_sd == 0], collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (perm.no < 99) {
    warning("At least 99 permutations are recommended for stable results.",
            call. = FALSE)
  }
  if (nrow(meta.dat) <= 40 && isTRUE(is.post.sample)) {
    if (verbose) {
      message("Posterior sampling is disabled when the sample size is <= 40.")
    }
    is.post.sample <- FALSE
  }
  if (!is.post.sample) {
    depth_model <- stats::lm(
      sqrt(colSums(feature.dat)) ~ meta.dat[[grp.name]]
    )
    depth_p <- stats::coef(summary(depth_model))[-1, "Pr(>|t|)"]
    if (length(depth_p) && min(depth_p, na.rm = TRUE) < 0.05) {
      warning(
        "Sequencing depth is associated with the variable of interest; ",
        "consider whether rarefaction or additional sensitivity analysis is needed.",
        call. = FALSE
      )
    }
  }

  relative <- sweep(feature.dat, 2L, colSums(feature.dat), "/")
  filter_ind <- rowMeans(relative != 0) >= prev.filter &
    rowMeans(relative) >= mean.abund.filter &
    matrixStats::rowMaxs(relative) >= max.abund.filter
  names(filter_ind) <- rownames(feature.dat)
  if (verbose) message(sum(!filter_ind), " features were filtered.")
  filter.features <- rownames(feature.dat)[!filter_ind]
  feature.dat <- feature.dat[filter_ind, , drop = FALSE]
  if (!nrow(feature.dat)) {
    stop("All features were removed by the ZicoSeq filters.", call. = FALSE)
  }

  keep_samples <- colSums(feature.dat) != 0
  if (any(!keep_samples)) {
    warning(
      sum(!keep_samples),
      " sample(s) were all zero after filtering and were removed.",
      call. = FALSE
    )
    feature.dat <- feature.dat[, keep_samples, drop = FALSE]
    meta.dat <- meta.dat[keep_samples, , drop = FALSE]
    raw_reference_total <- raw_reference_total[keep_samples]
    if (!is.null(strata)) strata <- strata[keep_samples]
  }

  sample_no <- ncol(feature.dat)
  otu_no <- nrow(feature.dat)
  feature_names <- rownames(feature.dat)
  depth <- colSums(feature.dat)
  if (verbose) {
    message("Testing ", otu_no, " features in ", sample_no, " samples.")
  }
  if (any(rowSums(feature.dat != 0) <= 2)) {
    warning(
      "Some features have fewer than three nonzero observations and little power.",
      call. = FALSE
    )
  }

  if (isTRUE(is.winsor)) {
    if (verbose) {
      message(
        "About ", round(length(depth) * outlier.pct),
        " outlier count(s) per feature will be replaced."
      )
    }
    winsorized <- vapply(seq_len(nrow(feature.dat)), function(i) {
      proportion <- feature.dat[i, ] / depth
      if (winsor.end == "top") {
        upper <- stats::quantile(proportion, 1 - outlier.pct)
        proportion[proportion >= upper] <- upper
      } else if (winsor.end == "bottom") {
        lower <- stats::quantile(proportion, outlier.pct)
        proportion[proportion <= lower] <- lower
      } else {
        upper <- stats::quantile(proportion, 1 - outlier.pct / 2)
        lower <- stats::quantile(proportion, outlier.pct / 2)
        proportion[proportion >= upper] <- upper
        proportion[proportion <= lower] <- lower
      }
      round(proportion * depth)
    }, numeric(ncol(feature.dat)))
    feature.dat <- t(winsorized)
    dimnames(feature.dat) <- list(feature_names, rownames(meta.dat))
    feature_sd <- matrixStats::rowSds(feature.dat)
    if (any(feature_sd == 0)) {
      stop(
        "Winsorization produced constant features: ",
        paste(feature_names[feature_sd == 0], collapse = ", "), ".",
        call. = FALSE
      )
    }
  }

  if (isTRUE(is.post.sample)) {
    if (verbose) message("Fitting beta mixtures.")
    posterior <- vapply(seq_len(nrow(feature.dat)), function(i) {
      count <- feature.dat[i, ]
      mixture <- tryCatch(
        .irs_bbmix_fit_mm(count, depth), error = function(e) NULL
      )
      if (!is.null(mixture)) {
        prop1 <- stats::rbeta(
          sample_no * post.sample.no,
          shape1 = count + mixture$shape1.1,
          shape2 = mixture$shape1.2 + depth - count
        )
        prop2 <- stats::rbeta(
          sample_no * post.sample.no,
          shape1 = count + mixture$shape2.1,
          shape2 = mixture$shape2.2 + depth - count
        )
        ifelse(
          stats::runif(sample_no * post.sample.no) <= mixture$q1,
          prop1, prop2
        )
      } else {
        proportion <- count / depth
        variance <- stats::var(proportion)
        average <- mean(proportion)
        a1 <- ((1 - average) / variance - 1 / average) * average^2
        a2 <- a1 * (1 / average - 1)
        if (!is.finite(a1) || a1 < 0) {
          stats::rbeta(
            sample_no * post.sample.no,
            shape1 = count + 1,
            shape2 = otu_no + depth - count
          )
        } else {
          stats::rbeta(
            sample_no * post.sample.no,
            shape1 = count + a1,
            shape2 = a2 + depth - count
          )
        }
      }
    }, numeric(sample_no * post.sample.no))
    posterior <- t(posterior)
    posterior_list <- lapply(seq_len(post.sample.no), function(i) {
      start <- (i - 1L) * sample_no + 1L
      posterior[, start:(start + sample_no - 1L), drop = FALSE]
    })
  } else {
    posterior_list <- list(sweep(feature.dat, 2L, depth, "/"))
    post.sample.no <- 1L
  }
  posterior_list <- lapply(posterior_list, function(x) {
    x[x <= min.prop] <- min.prop
    x
  })

  if (is.null(adj.name)) {
    M0 <- stats::model.matrix(~1, meta.dat)
  } else {
    adjustment_data <- meta.dat[, adj.name, drop = FALSE]
    if (anyNA(adjustment_data)) {
      stop("Adjustment variables cannot contain missing values.", call. = FALSE)
    }
    M0 <- stats::model.matrix(~., adjustment_data)
  }
  group_data <- meta.dat[, grp.name, drop = FALSE]
  if (anyNA(group_data)) {
    stop("The variable of interest cannot contain missing values.", call. = FALSE)
  }
  M1 <- stats::model.matrix(~., group_data)[, -1, drop = FALSE]
  if (!ncol(M1)) {
    stop("The variable of interest does not generate a testable coefficient.",
         call. = FALSE)
  }
  M01 <- cbind(M0, M1)
  qrX0 <- qr(M0, tol = 1e-7)
  Q0 <- qr.Q(qrX0)[, seq_len(qrX0$rank), drop = FALSE]
  H0 <- Q0 %*% t(Q0)
  qrX1 <- qr(M1, tol = 1e-7)
  Q1 <- qr.Q(qrX1)[, seq_len(qrX1$rank), drop = FALSE]
  qrX01 <- qr(M01, tol = 1e-7)
  Q01 <- qr.Q(qrX01)[, seq_len(qrX01$rank), drop = FALSE]
  R0 <- as.matrix(stats::resid(stats::lm(Q1 ~ Q0 - 1)))
  df.model <- ncol(Q01) - ncol(Q0)
  df.residual <- sample_no - ncol(Q01)
  if (df.model <= 0 || df.residual <= 0) {
    stop("The model has insufficient residual degrees of freedom.", call. = FALSE)
  }

  divisor <- raw_reference_total / depth
  zero_samples <- names(divisor)[!is.finite(divisor) | divisor <= 0]
  if (length(zero_samples)) {
    if (zero_ref_action == "error") {
      stop(
        "The original IRS reference total is zero for sample(s): ",
        paste(zero_samples, collapse = ", "), ".", call. = FALSE
      )
    }
    positive <- divisor[is.finite(divisor) & divisor > 0]
    if (!length(positive)) {
      stop("All IRS reference divisors are zero.", call. = FALSE)
    }
    divisor[!is.finite(divisor) | divisor <= 0] <- 0.5 * min(positive)
  }
  normalization_total <- divisor * depth

  if (verbose) message("Permutation testing.")
  func_no <- length(link.func)
  Y <- matrix(NA_real_, sample_no, func_no * otu_no * post.sample.no)
  for (k in seq_len(post.sample.no)) {
    for (j in seq_len(func_no)) {
      columns <- (k - 1L) * func_no * otu_no +
        func_no * (0:(otu_no - 1L)) + j
      Y[, columns] <- link.func[[j]](t(posterior_list[[k]]) / divisor)
    }
  }
  Y <- t(Y)
  total_ss <- rowSums(Y^2)
  model_ss <- rowSums((Y %*% Q01)^2)
  null_ss <- rowSums((Y %*% Q0)^2)
  explained_ss <- model_ss - null_ss
  residual_ss <- total_ss - model_ss

  get_permute_matrix <- utils::getFromNamespace("getPermuteMatrix", "vegan")
  permutation_index <- get_permute_matrix(
    perm.no, sample_no, strata = strata
  )
  perm.no <- nrow(permutation_index)
  permuted_ss <- vapply(seq_len(perm.no), function(ii) {
    if (verbose && ii %% 10L == 0L) cat(".")
    permuted_R <- R0[permutation_index[ii, ], , drop = FALSE]
    permuted_R <- permuted_R - H0 %*% permuted_R
    qr_permuted <- qr(permuted_R, tol = 1e-7)
    Q1_permuted <- qr.Q(qr_permuted)[, seq_len(qr_permuted$rank), drop = FALSE]
    permuted_model_ss <- null_ss + rowSums((Y %*% Q1_permuted)^2)
    c(permuted_model_ss - null_ss, total_ss - permuted_model_ss)
  }, numeric(2L * length(total_ss)))
  if (verbose) cat("\n")

  unit <- func_no * otu_no * post.sample.no
  permuted_explained <- permuted_ss[seq_len(unit), , drop = FALSE]
  permuted_residual <- permuted_ss[unit + seq_len(unit), , drop = FALSE]
  RSS.m <- array(residual_ss, c(func_no, otu_no, post.sample.no))
  RSS.m <- t(apply(RSS.m, c(1, 2), mean))
  F0.m <- array(
    (explained_ss / df.model) / (residual_ss / df.residual),
    c(func_no, otu_no, post.sample.no)
  )
  F0.m <- t(apply(F0.m, c(1, 2), mean))
  R2.m <- array(explained_ss / total_ss, c(func_no, otu_no, post.sample.no))
  R2.m <- t(apply(R2.m, c(1, 2), mean))
  F0 <- array(
    (explained_ss / df.model) / (residual_ss / df.residual),
    c(func_no, otu_no, post.sample.no)
  )
  Fp <- array(
    (permuted_explained / df.model) / (permuted_residual / df.residual),
    c(func_no, otu_no, post.sample.no, perm.no)
  )
  F0 <- apply(F0, c(1, 2), mean)
  Fp <- apply(Fp, c(1, 2, 4), mean)
  F0 <- apply(F0, 2, stats.combine.func)
  Fp <- apply(Fp, c(2, 3), stats.combine.func)

  if (mean(is.na(F0)) >= 0.1) {
    warning("More than 10% of observed F statistics are missing.", call. = FALSE)
  }
  if (mean(is.na(Fp)) >= 0.1) {
    warning("More than 10% of permuted F statistics are missing.", call. = FALSE)
  }
  na_ind <- is.na(F0)
  observed_complete <- F0[!na_ind]
  permuted_complete <- Fp[!na_ind, , drop = FALSE]
  p_raw_complete <- rowMeans(
    cbind(permuted_complete, observed_complete) >= observed_complete
  )
  p_fdr_complete <- .irs_perm_fdr_adj(
    observed_complete, permuted_complete
  )
  p.raw <- .irs_na_pad(p_raw_complete, na_ind)
  p.adj.fdr <- .irs_na_pad(p_fdr_complete, na_ind)
  names(p.raw) <- names(p.adj.fdr) <- feature_names
  rownames(R2.m) <- rownames(RSS.m) <- rownames(F0.m) <- feature_names
  colnames(R2.m) <- colnames(RSS.m) <- colnames(F0.m) <- paste0(
    "Func", seq_len(func_no)
  )
  if (isTRUE(is.fwer)) {
    p_fwer_complete <- .irs_perm_fwer_adj(
      observed_complete, permuted_complete
    )
    p.adj.fwer <- .irs_na_pad(p_fwer_complete, na_ind)
    names(p.adj.fwer) <- feature_names
  } else {
    p.adj.fwer <- NULL
  }

  normalized_feature_data <- sweep(
    feature.dat, 2L, normalization_total, "/"
  )
  coefficient_list <- lapply(link.func, function(fun) {
    solve(crossprod(M01), crossprod(M01, t(fun(normalized_feature_data))))
  })
  returned_feature_data <- if (isTRUE(return.feature.dat)) feature.dat else NULL
  if (verbose) message("Completed.")

  list(
    call = this_call,
    feature.dat = returned_feature_data,
    meta.dat = meta.dat,
    grp.name = grp.name,
    filter.features = filter.features,
    ref.features = ref.features,
    R2 = R2.m,
    F0 = F0.m,
    RSS = RSS.m,
    df.model = df.model,
    df.residual = df.residual,
    coef.list = coefficient_list,
    p.raw = p.raw,
    p.adj.fdr = p.adj.fdr,
    p.adj.fwer = p.adj.fwer,
    irs = list(
      integration_mode = "normalization",
      reference_taxa = ref.features,
      raw_reference_total = raw_reference_total,
      normalization_total = normalization_total,
      zero_reference_samples = zero_samples,
      denominator_source = "original_counts",
      zicoseq_source = "GUniFrac 1.9"
    )
  )
}
