#' Transform Mixed Data to Uniform Pseudo-Observations
#'
#' Converts a \code{structureR_data} object (or a plain numeric
#' \code{data.frame}/\code{matrix}) to the \eqn{(0,1)} hypercube by applying
#' an Empirical Cumulative Distribution Function (ECDF) / Probability
#' Integral Transform (PIT) to each column. Gene expression columns that
#' exhibit zero-inflation are handled with a two-component mixture approach:
#' zeros receive mass on \eqn{[0, p_0]} and positive values are mapped via
#' their conditional ECDF onto \eqn{(p_0, 1)}, where \eqn{p_0} is the
#' empirical zero-proportion.
#'
#' @param vine_data Either a \code{"structureR_data"} object returned by
#'   \code{\link{PrepVineData}}, or a numeric \code{data.frame} /
#'   \code{matrix} whose columns are to be transformed.
#' @param zero_augmented_genes Optional character vector of column names that
#'   should be transformed with the zero-augmented PIT. Only meaningful when
#'   \code{vine_data} is a \code{data.frame}/\code{matrix}. When
#'   \code{vine_data} is a \code{structureR_data} object, gene expression
#'   columns are automatically treated as zero-augmented.
#' @param use_rcpp Logical. Use the Rcpp-accelerated ECDF implementation for
#'   non-zero-augmented columns (default: \code{TRUE}).
#'
#' @return A numeric matrix of pseudo-observations in \eqn{(0,1)}, one row
#'   per observation and one column per variable. The matrix preserves column
#'   names.
#'
#' @examples
#' set.seed(1)
#' df <- data.frame(x = rnorm(50), y = rbeta(50, 2, 5), z = rpois(50, 1))
#' pseudo <- TransformMarginals(df)
#' stopifnot(all(pseudo > 0 & pseudo < 1))
#'
#' @export
TransformMarginals <- function(vine_data,
                               zero_augmented_genes = NULL,
                               use_rcpp             = TRUE) {

  ## ---- 1. Extract a numeric matrix from the input -------------------------
  if (inherits(vine_data, "structureR_data")) {
    # Merge proportions and subject metadata on the unit key
    unit_var <- vine_data$params$unit_var
    merged   <- vine_data$proportions

    if (!is.null(vine_data$subject_metadata)) {
      merged <- merge(merged, vine_data$subject_metadata,
                      by = unit_var, all.x = TRUE)
    }

    if (!is.null(vine_data$gene_expression)) {
      merged <- merge(merged, vine_data$gene_expression,
                      by = unit_var, all.x = TRUE)
      # Mark gene columns for zero-augmented treatment
      gene_cols            <- setdiff(colnames(vine_data$gene_expression), unit_var)
      zero_augmented_genes <- union(zero_augmented_genes, gene_cols)
    }

    # Drop the unit-ID column and any non-numeric columns
    merged[[unit_var]] <- NULL
    num_flags    <- vapply(merged, is.numeric, logical(1))
    data_matrix  <- merged[, num_flags, drop = FALSE]

  } else {
    # Plain data.frame or matrix
    if (is.matrix(vine_data)) {
      data_matrix <- as.data.frame(vine_data)
    } else {
      data_matrix <- vine_data
    }
    num_flags   <- vapply(data_matrix, is.numeric, logical(1))
    data_matrix <- data_matrix[, num_flags, drop = FALSE]
  }

  if (ncol(data_matrix) == 0L) {
    stop("No numeric columns found in 'vine_data'.")
  }

  ## ---- 2. Transform each column -------------------------------------------
  n          <- nrow(data_matrix)
  pseudo_obs <- data_matrix   # copy; will overwrite values

  for (col_nm in colnames(data_matrix)) {
    x <- data_matrix[[col_nm]]

    if (!is.null(zero_augmented_genes) && col_nm %in% zero_augmented_genes) {
      pseudo_obs[[col_nm]] <- .zero_augmented_pit(x)
    } else if (use_rcpp) {
      pseudo_obs[[col_nm]] <- rcpp_ecdf_transform(x)
    } else {
      # Pure R rank-based pseudo-observations
      pseudo_obs[[col_nm]] <- rank(x, ties.method = "average") / (n + 1.0)
    }
  }

  as.matrix(pseudo_obs)
}


## ---------------------------------------------------------------------------
## Internal helpers
## ---------------------------------------------------------------------------

#' Zero-Augmented Probability Integral Transform
#'
#' @param x Numeric vector (possibly zero-inflated).
#' @return Numeric vector in (0, 1).
#' @noRd
.zero_augmented_pit <- function(x) {
  n          <- length(x)
  prop_zero  <- mean(x == 0)
  u          <- numeric(n)
  is_zero    <- (x == 0)

  if (prop_zero == 0) {
    # No zeros — standard ECDF
    return(rank(x, ties.method = "average") / (n + 1.0))
  }

  # Zeros: random uniform draw on [epsilon, prop_zero]
  eps         <- .Machine$double.eps
  u[is_zero]  <- stats::runif(sum(is_zero),
                               min = eps,
                               max = max(prop_zero, 1e-10))

  # Non-zeros: conditional ECDF mapped onto (prop_zero, 1)
  if (any(!is_zero)) {
    nz_vals       <- x[!is_zero]
    ecdf_nz       <- stats::ecdf(nz_vals)
    nz_pseudo     <- ecdf_nz(nz_vals)          # in [0, 1]
    # Shift into (prop_zero, 1)
    u[!is_zero]   <- prop_zero + (1.0 - prop_zero) * nz_pseudo
  }

  # Clamp strictly inside (0, 1)
  pmax(pmin(u, 1.0 - .Machine$double.eps), .Machine$double.eps)
}
