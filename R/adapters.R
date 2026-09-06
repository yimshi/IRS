.require_optional_package <- function(package, purpose) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(
      "Package '", package, "' is required for ", purpose,
      ". Install it and try again.", call. = FALSE
    )
  }
}

.aligned_fit_vector <- function(fit, sample_names, value, object_name) {
  .check_irs_fit(fit)
  if (is.null(sample_names) || anyDuplicated(sample_names)) {
    stop(object_name, " must have unique sample names.", call. = FALSE)
  }
  missing <- setdiff(sample_names, fit$sample_names)
  extra <- setdiff(fit$sample_names, sample_names)
  if (length(missing) || length(extra)) {
    stop(
      object_name, " samples do not match the IRS fit. Missing from fit: ",
      paste(missing, collapse = ", "), "; missing from object: ",
      paste(extra, collapse = ", "), ".", call. = FALSE
    )
  }
  value[sample_names]
}

#' Inject IRS size factors into a DESeq2 data set
#'
#' The original integer counts are retained. The injected size factors are the
#' geometric-mean-centered IRS reference totals, so DESeq2 normalized counts
#' are a common constant multiple of the IRS ratios.
#'
#' @param dds A `DESeqDataSet`.
#' @param fit An `irs_fit` object.
#' @param overwrite Replace existing DESeq2 size factors. Normalization-factor
#'   matrices are never silently removed.
#'
#' @return The modified `DESeqDataSet`.
#' @export
inject_deseq2 <- function(dds, fit, overwrite = FALSE) {
  .require_optional_package("DESeq2", "inject_deseq2()")
  if (!inherits(dds, "DESeqDataSet")) {
    stop("dds must inherit from DESeqDataSet.", call. = FALSE)
  }
  nf <- DESeq2::normalizationFactors(dds)
  if (!is.null(nf)) {
    stop(
      "dds already contains a normalization-factor matrix, which overrides ",
      "size factors. Remove it explicitly before injecting IRS.", call. = FALSE
    )
  }
  existing <- DESeq2::sizeFactors(dds)
  if (!is.null(existing) && !isTRUE(overwrite)) {
    stop("dds already has size factors; use overwrite = TRUE to replace them.",
         call. = FALSE)
  }
  sample_names <- colnames(DESeq2::counts(dds))
  sf <- .aligned_fit_vector(fit, sample_names, fit$size_factors, "dds")
  DESeq2::sizeFactors(dds) <- unname(sf)
  dds
}

#' Inject IRS effective library sizes into an edgeR DGEList
#'
#' Sets edgeR normalization factors so each effective library size is a common
#' multiple of the IRS reference total. Original integer counts and native
#' edgeR model fitting are retained.
#'
#' @param y An edgeR `DGEList`.
#' @param fit An `irs_fit` object.
#' @param overwrite Replace non-unity normalization factors.
#'
#' @return The modified `DGEList`.
#' @export
inject_edger <- function(y, fit, overwrite = FALSE) {
  .require_optional_package("edgeR", "inject_edger()")
  if (!inherits(y, "DGEList")) {
    stop("y must inherit from DGEList.", call. = FALSE)
  }
  if (!is.null(y$offset)) {
    stop(
      "y already contains an offset matrix, which would override effective ",
      "library sizes. Remove it explicitly before injecting IRS.", call. = FALSE
    )
  }
  current <- y$samples$norm.factors
  if (!is.null(current) && any(abs(current - 1) > 1e-12) && !isTRUE(overwrite)) {
    stop(
      "y already has non-unity normalization factors; use overwrite = TRUE ",
      "to replace them.", call. = FALSE
    )
  }
  sample_names <- colnames(y$counts)
  ref_total <- .aligned_fit_vector(
    fit, sample_names, fit$reference_total, "y"
  )
  lib_size <- stats::setNames(y$samples$lib.size, sample_names)
  if (any(!is.finite(lib_size)) || any(lib_size <= 0)) {
    stop("edgeR library sizes must be finite and positive.", call. = FALSE)
  }
  scaling_constant <- .geometric_mean(lib_size) / .geometric_mean(ref_total)
  norm_factors <- scaling_constant * ref_total / lib_size
  y$samples$norm.factors <- unname(norm_factors)
  y
}

#' Prepare IRS-normalized input for MaAsLin2
#'
#' MaAsLin2 should be called with `normalization = "NONE"` and, for the
#' currently validated integration, `analysis_method = "LM"`. Because current
#' MaAsLin2 releases filter before their internal normalization step, non-zero
#' abundance or variance thresholds can have different semantics when external
#' normalized input is supplied.
#'
#' @param fit An `irs_fit` object.
#' @param counts Optional taxa-by-samples count matrix.
#'
#' @return A list containing `input_data` and recommended native arguments.
#' @export
prepare_maaslin2 <- function(fit, counts = NULL) {
  input <- normalized_counts(
    fit, counts, orientation = "samples_by_taxa"
  )
  out <- list(
    input_data = as.data.frame(input, check.names = FALSE),
    normalization = "NONE",
    analysis_method = "LM",
    integration_mode = "normalization",
    validated_scope = "MaAsLin2 linear model"
  )
  class(out) <- c("irs_maaslin2_input", "list")
  out
}

#' Return IRS reference indices for DACOMP
#'
#' @param fit An `irs_fit` object.
#' @param X DACOMP count matrix with samples in rows and taxa in columns.
#'
#' @return Integer reference-taxon column indices. This is a reference-set
#'   integration; DACOMP retains its own normalization and testing definition.
#' @export
reference_for_dacomp <- function(fit, X) {
  X <- as.matrix(X)
  if (is.null(colnames(X))) {
    stop("X must have taxon column names.", call. = FALSE)
  }
  unname(reference_indices(fit, X, taxa_are_rows = FALSE))
}

#' Return IRS denominator indices for ALDEx2
#'
#' ALDEx2 removes all-zero features before interpreting numeric denominator
#' indices. This helper performs the same filtering before matching IRS taxa.
#'
#' @param fit An `irs_fit` object.
#' @param reads ALDEx2 count matrix with taxa in rows and samples in columns.
#'
#' @return Integer indices for the `denom` argument to [ALDEx2::aldex.clr()].
#'   This is a reference-set integration: ALDEx2 uses its own geometric-mean
#'   denominator rather than the IRS reference sum.
#' @export
denom_for_aldex2 <- function(fit, reads) {
  reads <- .validate_irs_counts(
    reads, require_names = TRUE, integerish = TRUE, object_name = "reads"
  )
  kept <- rowSums(reads) > 0
  kept_names <- rownames(reads)[kept]
  missing <- setdiff(fit$reference_taxa, kept_names)
  if (length(missing)) {
    stop(
      "IRS reference taxa would be removed by ALDEx2 as all-zero features: ",
      paste(missing, collapse = ", "), ".", call. = FALSE
    )
  }
  unname(match(fit$reference_taxa, kept_names))
}

#' Prepare IRS-normalized input for the frequency-scale LDM
#'
#' @param fit An `irs_fit` object.
#' @param counts Optional taxa-by-samples count matrix.
#'
#' @return A list containing a samples-by-taxa OTU table and the native LDM
#'   arguments required to avoid a second normalization step.
#' @export
prepare_ldm <- function(fit, counts = NULL) {
  otu <- normalized_counts(
    fit, counts, orientation = "samples_by_taxa"
  )
  out <- list(
    otu_table = as.data.frame(otu, check.names = FALSE),
    scale.otu.table = FALSE,
    freq.scale.only = TRUE,
    comp.anal = FALSE,
    integration_mode = "normalization",
    validated_scope = "LDM frequency-scale analysis"
  )
  class(out) <- c("irs_ldm_input", "list")
  out
}

#' Taxon-wise Wilcoxon tests after IRS normalization
#'
#' @param fit An `irs_fit` object.
#' @param counts Optional taxa-by-samples count matrix.
#' @param group Two-level grouping variable aligned to the fitted samples, or a
#'   metadata column name when `metadata` is supplied.
#' @param metadata Optional sample metadata with matching row names.
#' @param p_adjust Multiple-testing method passed to [stats::p.adjust()].
#' @param exact Passed to [stats::wilcox.test()].
#'
#' @return A data frame with taxon-wise statistics and adjusted p-values.
#' @export
irs_wilcox <- function(fit, counts = NULL, group, metadata = NULL,
                       p_adjust = "BH", exact = FALSE) {
  x <- normalized_counts(fit, counts, orientation = "taxa_by_samples")
  if (length(group) == 1L && is.character(group) && !is.null(metadata)) {
    metadata <- .align_irs_metadata(x, metadata)
    if (!group %in% colnames(metadata)) {
      stop("group was not found in metadata.", call. = FALSE)
    }
    group <- metadata[[group]]
  }
  if (length(group) != ncol(x) || anyNA(group)) {
    stop("group must contain one non-missing value per sample.", call. = FALSE)
  }
  group <- factor(group)
  if (nlevels(group) != 2L) {
    stop("irs_wilcox() currently requires exactly two groups.", call. = FALSE)
  }

  tests <- lapply(seq_len(nrow(x)), function(i) {
    values <- x[i, ]
    if (length(unique(values)) < 2L) {
      return(c(statistic = 0, p_value = 1))
    }
    ans <- suppressWarnings(stats::wilcox.test(values ~ group, exact = exact))
    c(statistic = unname(ans$statistic), p_value = ans$p.value)
  })
  tests <- do.call(rbind, tests)
  data.frame(
    taxon = rownames(x),
    statistic = tests[, "statistic"],
    p_value = tests[, "p_value"],
    p_adjusted = stats::p.adjust(tests[, "p_value"], method = p_adjust),
    integration_mode = "normalization",
    row.names = rownames(x),
    check.names = FALSE
  )
}
