#' Vandeputte 2017 Crohn's disease microbiome data
#'
#' Genus-level 16S rRNA gene count data and flow-cytometry microbial-load
#' measurements for 95 samples from a Crohn's disease case-control cohort.
#' The object is provided for examples of IRS normalization and external
#' comparison with quantitative microbiome profiling.
#'
#' @format A list with two elements:
#' \describe{
#'   \item{counts}{An integer matrix with 166 genera in rows and 95 samples in
#'   columns.}
#'   \item{sample_data}{A data frame with one row per sample and the variables
#'   `sample_id`, `ena_id`, `group`, and `microbial_load`. The load is the
#'   average flow-cytometry cell count per gram of frozen feces.}
#' }
#' @source
#' Vandeputte D, Kathagen G, D'hoe K, et al. (2017). Quantitative microbiome
#' profiling links gut community variation to microbial load. *Nature* 551,
#' 507--511. \doi{10.1038/nature24460}. Processed count and load files were
#' obtained from the CellCount_Nishijima_2024 data release
#' (\doi{10.5281/zenodo.14243685}) and are redistributed under
#' \href{https://creativecommons.org/licenses/by/4.0/}{CC BY 4.0}.
#' @examples
#' data(vandeputte2017)
#' dim(vandeputte2017$counts)
#' table(vandeputte2017$sample_data$group)
"vandeputte2017"

#' Paired DAA benchmark summaries
#'
#' Precomputed summaries comparing native and IRS-integrated normalization for
#' eight differential abundance analysis workflows. Results are included for
#' the simulation used in the getting-started vignette and for the Vandeputte
#' 2017 Crohn's disease data.
#'
#' @format A data frame with 16 rows and the following variables:
#' \describe{
#'   \item{dataset}{Benchmark data set.}
#'   \item{method}{Downstream DAA workflow.}
#'   \item{native_variant}{Normalization used in the native comparison.}
#'   \item{native_discoveries, irs_discoveries}{Numbers of adjusted-p-value or
#'   q-value discoveries at 0.05.}
#'   \item{native_benchmark_positive, irs_benchmark_positive}{Discoveries that
#'   match simulated signals or the Vandeputte QMP-positive set.}
#'   \item{native_benchmark_negative, irs_benchmark_negative}{Discoveries
#'   outside the corresponding benchmark-positive set.}
#'   \item{native_recall, irs_recall, native_f1, irs_f1}{Concordance metrics
#'   relative to the benchmark-positive set.}
#'   \item{analysis_settings}{Package versions and principal Monte Carlo or
#'   permutation settings used for the benchmark.}
#' }
#' @details
#' The simulation has five known differential taxa. For the Vandeputte data,
#' benchmark-positive taxa are significant in a Wilcoxon comparison of
#' flow-cytometry-derived quantitative microbiome profiles after BH adjustment
#' at 0.05.
#' @examples
#' data(irs_daa_benchmark)
#' subset(irs_daa_benchmark, dataset == "simulation")
"irs_daa_benchmark"
