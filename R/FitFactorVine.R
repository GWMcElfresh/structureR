#' Fit a Factor-Vine Copula Model
#'
#' Fits a vine copula using \pkg{rvinecopulib} with a "factor-first"
#' structure: the dominant covariate (\code{factor_var}) is placed as the
#' root of a D-vine so that all pair-copulas involving it are estimated
#' first. The remaining dependency structure is discovered automatically by
#' \code{\link[rvinecopulib]{vinecop}}.
#'
#' @param pseudo_obs A numeric matrix of pseudo-observations in \eqn{(0,1)},
#'   typically the output of \code{\link{TransformMarginals}}.
#' @param factor_var Either a character string matching one of
#'   \code{colnames(pseudo_obs)}, or a positive integer giving the column
#'   index of the factor variable to use as the root. If \code{NULL}, column
#'   1 is used.
#' @param var_names Optional character vector of variable names. If provided,
#'   it overrides \code{colnames(pseudo_obs)}.
#' @param ... Additional arguments forwarded to
#'   \code{\link[rvinecopulib]{vinecop}}, e.g. \code{family_set} or
#'   \code{selcrit}.
#'
#' @return An S3 object of class \code{"structureR_vine"} with components:
#'   \describe{
#'     \item{vine}{The fitted \code{vinecop} object from
#'       \pkg{rvinecopulib}.}
#'     \item{marginals}{A list describing the marginal setup: pseudo-obs,
#'       column ordering, and variable names.}
#'     \item{tree_plot}{A \pkg{ggplot2} plot of the vine tree structure, or
#'       \code{NULL} if \pkg{ggplot2} is not available.}
#'     \item{params}{A list of model-level metadata.}
#'   }
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' pseudo <- matrix(runif(200), nrow = 50, ncol = 4)
#' colnames(pseudo) <- c("disease_stage", "prop_0", "prop_1", "age")
#' fit <- FitFactorVine(pseudo, factor_var = "disease_stage")
#' print(fit)
#' }
#'
#' @export
FitFactorVine <- function(pseudo_obs,
                          factor_var = NULL,
                          var_names  = NULL,
                          ...) {

  if (!requireNamespace("rvinecopulib", quietly = TRUE)) {
    stop("The 'rvinecopulib' package is required. Install it with: ",
         "install.packages('rvinecopulib')")
  }

  if (!is.matrix(pseudo_obs)) {
    pseudo_obs <- as.matrix(pseudo_obs)
  }
  if (any(pseudo_obs <= 0) || any(pseudo_obs >= 1)) {
    stop("'pseudo_obs' must contain values strictly in (0, 1). ",
         "Run TransformMarginals() first.")
  }

  n_vars <- ncol(pseudo_obs)
  if (!is.null(var_names)) {
    colnames(pseudo_obs) <- var_names
  }

  ## ---- Determine factor (root) column index --------------------------------
  if (!is.null(factor_var) && is.character(factor_var)) {
    factor_idx <- which(colnames(pseudo_obs) == factor_var)
    if (length(factor_idx) == 0L) {
      warning("'factor_var' = '", factor_var, "' not found in column names. ",
              "Using column 1 as the factor.")
      factor_idx <- 1L
    } else {
      factor_idx <- factor_idx[[1L]]
    }
  } else if (is.numeric(factor_var) && length(factor_var) == 1L) {
    factor_idx <- as.integer(factor_var)
  } else {
    factor_idx <- 1L
  }

  ## ---- Re-order columns: factor first, then the rest ----------------------
  other_idx      <- setdiff(seq_len(n_vars), factor_idx)
  ordered_idx    <- c(factor_idx, other_idx)
  pseudo_ordered <- pseudo_obs[, ordered_idx, drop = FALSE]

  ## ---- Fit D-vine with factor as root (first node) ------------------------
  fitted_vine <- tryCatch({
    structure_dvine <- rvinecopulib::dvine_structure(order = seq_len(n_vars))
    rvinecopulib::vinecop(data      = pseudo_ordered,
                          structure = structure_dvine,
                          ...)
  }, error = function(e) {
    warning("D-vine fitting failed (", conditionMessage(e), "); ",
            "falling back to automatic R-vine structure selection.")
    rvinecopulib::vinecop(data = pseudo_ordered, ...)
  })

  ## ---- Optional dependency-tree visualisation ------------------------------
  tree_plot <- NULL
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    tree_plot <- tryCatch(plot(fitted_vine), error = function(e) NULL)
  }

  ## ---- Assemble S3 result --------------------------------------------------
  structure(
    list(
      vine = fitted_vine,
      marginals = list(
        pseudo_obs     = pseudo_obs,
        ordered_pseudo = pseudo_ordered,
        ordered_idx    = ordered_idx,
        factor_var     = if (is.null(factor_var)) colnames(pseudo_obs)[factor_idx]
                         else factor_var,
        factor_idx     = factor_idx,
        var_names      = colnames(pseudo_obs)
      ),
      tree_plot = tree_plot,
      params    = list(
        n_vars     = n_vars,
        n_obs      = nrow(pseudo_obs),
        factor_var = if (is.null(factor_var)) colnames(pseudo_obs)[factor_idx]
                     else factor_var
      )
    ),
    class = "structureR_vine"
  )
}


## ---------------------------------------------------------------------------
## S3 methods
## ---------------------------------------------------------------------------

#' @export
print.structureR_vine <- function(x, ...) {
  cat("structureR Factor-Vine Model\n")
  cat("-----------------------------\n")
  cat("Variables   :", x$params$n_vars, "\n")
  cat("Observations:", x$params$n_obs,  "\n")
  if (!is.null(x$params$factor_var)) {
    cat("Factor root :", x$params$factor_var, "\n")
  }
  if (!is.null(x$marginals$var_names)) {
    cat("Var names   :", paste(x$marginals$var_names, collapse = ", "), "\n")
  }
  invisible(x)
}

#' @export
summary.structureR_vine <- function(object, ...) {
  print(object)
  cat("\nVine structure:\n")
  print(object$vine)
  invisible(object)
}
