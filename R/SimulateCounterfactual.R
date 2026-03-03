#' Generate Counterfactual Samples from a Fitted Factor-Vine
#'
#' Fixes the factor (disease) node at a user-specified pseudo-observation
#' value and samples from the remaining conditional distribution using the
#' inverse Rosenblatt transform provided by \pkg{rvinecopulib}. The result
#' represents synthetic cell proportions and gene expression profiles that
#' would be expected under the specified disease state.
#'
#' @param fitted_vine A \code{"structureR_vine"} object returned by
#'   \code{\link{FitFactorVine}}.
#' @param factor_value A numeric scalar in \eqn{(0, 1)}. The value at which
#'   the factor (disease) variable is fixed on the pseudo-observation
#'   (uniform) scale.
#' @param n_sim Positive integer. Number of synthetic observations to
#'   generate (default: 100).
#'
#' @return An S3 object of class \code{"structureR_counterfactual"} with
#'   components:
#'   \describe{
#'     \item{samples}{A numeric matrix (\code{n_sim} x \code{n_vars}) of
#'       simulated pseudo-observations. Column 1 is the fixed factor
#'       variable; remaining columns are conditionally simulated.}
#'     \item{factor_var}{Name of the factor variable (may be \code{NULL}).}
#'     \item{factor_value}{The fixed value used for the factor variable.}
#'     \item{n_sim}{Number of simulations requested.}
#'   }
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' pseudo <- matrix(runif(200), nrow = 50, ncol = 4)
#' colnames(pseudo) <- c("disease_stage", "prop_0", "prop_1", "age")
#' fit <- FitFactorVine(pseudo, factor_var = "disease_stage")
#' cf  <- SimulateCounterfactual(fit, factor_value = 0.9, n_sim = 50)
#' print(cf)
#' }
#'
#' @export
SimulateCounterfactual <- function(fitted_vine,
                                   factor_value,
                                   n_sim = 100L) {

  if (!inherits(fitted_vine, "structureR_vine")) {
    stop("'fitted_vine' must be a 'structureR_vine' object from FitFactorVine().")
  }
  if (!requireNamespace("rvinecopulib", quietly = TRUE)) {
    stop("The 'rvinecopulib' package is required. Install it with: ",
         "install.packages('rvinecopulib')")
  }

  if (!is.numeric(factor_value) || length(factor_value) != 1L ||
      factor_value <= 0 || factor_value >= 1) {
    stop("'factor_value' must be a single numeric value strictly in (0, 1).")
  }

  n_sim  <- as.integer(n_sim)
  n_vars <- fitted_vine$params$n_vars

  ## ---- Conditional simulation via inverse Rosenblatt ----------------------
  # Build a matrix where column 1 (factor) is fixed and the remaining
  # columns are drawn uniformly: [factor_value, U_2, ..., U_p].
  # Applying the inverse Rosenblatt transform converts the uniforms to
  # samples from the conditional copula C(u_2,...,u_p | u_1 = factor_value).

  u_input        <- matrix(stats::runif(n_sim * n_vars),
                            nrow = n_sim, ncol = n_vars)
  u_input[, 1L]  <- factor_value   # fix the factor column

  simulated <- tryCatch(
    rvinecopulib::inverse_rosenblatt(u_input, fitted_vine$vine),
    error = function(e) {
      warning("inverse_rosenblatt failed (", conditionMessage(e), "); ",
              "falling back to unconditional sampling with post-filtering.")
      # Fallback: oversample and pick rows closest to the target
      n_over   <- max(n_sim * 20L, 2000L)
      sim_full <- rvinecopulib::rvinecop(n_over, fitted_vine$vine)
      dist_vec <- abs(sim_full[, 1L] - factor_value)
      top_idx  <- order(dist_vec)[seq_len(min(n_sim, nrow(sim_full)))]
      sim_full[top_idx, , drop = FALSE]
    }
  )

  ## ---- Column names -------------------------------------------------------
  ordered_names <- fitted_vine$marginals$var_names[fitted_vine$marginals$ordered_idx]
  if (!is.null(ordered_names) &&
      length(ordered_names) == ncol(simulated)) {
    colnames(simulated) <- ordered_names
  }

  ## ---- Result -------------------------------------------------------------
  structure(
    list(
      samples      = simulated,
      factor_var   = fitted_vine$marginals$factor_var,
      factor_value = factor_value,
      n_sim        = n_sim
    ),
    class = "structureR_counterfactual"
  )
}


## ---------------------------------------------------------------------------
## S3 methods
## ---------------------------------------------------------------------------

#' @export
print.structureR_counterfactual <- function(x, ...) {
  cat("structureR Counterfactual Simulation\n")
  cat("--------------------------------------\n")
  if (!is.null(x$factor_var)) {
    cat("Factor     :", x$factor_var, "\n")
  }
  cat("Fixed value:", x$factor_value, "\n")
  cat("Simulations:", x$n_sim, "\n")
  cat("Variables  :", ncol(x$samples), "\n")
  invisible(x)
}
