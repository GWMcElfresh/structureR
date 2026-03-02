#' structureR: Factor-Vine Hierarchical Copula Modeling for Single-Cell Data
#'
#' Models multi-level dependencies in single-cell data using Factor-Vine
#' Hierarchical Copulas. The package extracts three distinct data layers from
#' a Seurat Object — subject-level metadata, computed cellular proportions,
#' and gene expression — and fits a hierarchical dependency structure that
#' isolates a dominant covariate (e.g., disease state) as the root factor.
#'
#' Main functions:
#' \describe{
#'   \item{\code{\link{PrepVineData}}}{Extract and structure data from a
#'     Seurat object or metadata data.frame.}
#'   \item{\code{\link{TransformMarginals}}}{Transform mixed data types to
#'     uniform pseudo-observations using ECDF/PIT.}
#'   \item{\code{\link{FitFactorVine}}}{Fit a factor-vine copula model with
#'     the dominant covariate as the root node.}
#'   \item{\code{\link{SimulateCounterfactual}}}{Generate synthetic samples
#'     by fixing the disease node at a target value.}
#' }
#'
#' @docType package
#' @name structureR-package
#' @aliases structureR
#' @useDynLib structureR, .registration = TRUE
#' @importFrom Rcpp sourceCpp
"_PACKAGE"
