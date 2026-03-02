#' Prepare Data for Vine Copula Modeling
#'
#' Extracts and structures three data layers from a Seurat object (or a
#' metadata data.frame): subject-level metadata, computed cellular
#' proportions, and optionally gene expression. The resulting object is the
#' primary input to \code{\link{TransformMarginals}}.
#'
#' @param seuratObj A Seurat object. Either \code{seuratObj} or
#'   \code{metadata} must be provided.
#' @param metadata A \code{data.frame} containing per-cell metadata (e.g.,
#'   taken directly from \code{seuratObj@@meta.data}). Used when
#'   \code{seuratObj} is \code{NULL}.
#' @param cluster_var Character. Column name in metadata that holds the
#'   cluster assignment for each cell (default: \code{"seurat_clusters"}).
#' @param unit_var Character. Column name representing the unit of data
#'   collection, e.g. sequencing lane / cDNA library ID (default:
#'   \code{"cDNA_ID"}).
#' @param subject_var Character. Column name for the subject / donor
#'   identifier (default: \code{"subject_id"}).
#' @param disease_var Character. Column name for the primary disease /
#'   condition covariate (default: \code{"disease_stage"}).
#' @param gene_list Optional character vector of gene names to aggregate from
#'   the Seurat assay data. Only used when \code{seuratObj} is provided.
#'
#' @return An S3 object of class \code{"structureR_data"} with components:
#'   \describe{
#'     \item{proportions}{A \code{data.frame} of cellular proportions per
#'       unit (rows = units, columns = cluster proportions + \code{unit_var}).}
#'     \item{subject_metadata}{A \code{data.frame} of unique subject-level
#'       covariates, or \code{NULL} if \code{subject_var} is not found.}
#'     \item{gene_expression}{A \code{data.frame} of mean gene expression per
#'       unit (aggregated across cells), or \code{NULL} if not requested.}
#'     \item{params}{A named list of the parameter values used.}
#'   }
#'
#' @examples
#' set.seed(42)
#' meta <- data.frame(
#'   cDNA_ID        = rep(paste0("s", 1:4), each = 50),
#'   subject_id     = rep(paste0("subj", 1:4), each = 50),
#'   disease_stage  = rep(c("healthy", "disease", "healthy", "disease"), each = 50),
#'   seurat_clusters = sample(0:3, 200, replace = TRUE),
#'   age            = rep(c(25, 45, 30, 50), each = 50)
#' )
#' vd <- PrepVineData(metadata = meta)
#' print(vd)
#'
#' @export
PrepVineData <- function(seuratObj   = NULL,
                         metadata    = NULL,
                         cluster_var = "seurat_clusters",
                         unit_var    = "cDNA_ID",
                         subject_var = "subject_id",
                         disease_var = "disease_stage",
                         gene_list   = NULL) {

  ## ---- 1. Obtain per-cell metadata ----------------------------------------
  if (!is.null(seuratObj)) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("The 'Seurat' package is required when 'seuratObj' is provided. ",
           "Install it with: BiocManager::install('Seurat')")
    }
    meta <- seuratObj@meta.data
  } else if (!is.null(metadata)) {
    meta <- as.data.frame(metadata, stringsAsFactors = FALSE)
  } else {
    stop("Either 'seuratObj' or 'metadata' must be provided.")
  }

  ## ---- 2. Validate required columns ----------------------------------------
  required_cols <- c(cluster_var, unit_var)
  missing_cols  <- setdiff(required_cols, colnames(meta))
  if (length(missing_cols) > 0) {
    stop("The following required columns are missing from the metadata: ",
         paste(missing_cols, collapse = ", "))
  }

  ## ---- 3. Compute cellular proportions ------------------------------------
  # Count cells per (unit x cluster), then normalise per unit
  count_table   <- table(meta[[unit_var]], meta[[cluster_var]])
  prop_matrix   <- count_table / rowSums(count_table)
  proportions_df <- as.data.frame.matrix(prop_matrix, stringsAsFactors = FALSE)
  # Prefix cluster columns to avoid name clashes
  colnames(proportions_df) <- paste0("cluster_", colnames(proportions_df))
  proportions_df[[unit_var]] <- rownames(proportions_df)
  rownames(proportions_df)   <- NULL

  ## ---- 4. Extract subject-level metadata -----------------------------------
  subject_meta <- NULL
  if (subject_var %in% colnames(meta)) {
    keep_cols <- c(unit_var, subject_var)
    if (disease_var %in% colnames(meta)) {
      keep_cols <- c(keep_cols, disease_var)
    }
    # Append any additional numeric columns (e.g., age, BMI)
    other_cols   <- setdiff(colnames(meta), c(cluster_var, keep_cols))
    numeric_cols <- other_cols[
      vapply(meta[, other_cols, drop = FALSE], is.numeric, logical(1))
    ]
    keep_cols    <- unique(c(keep_cols, numeric_cols))
    subject_meta <- unique(meta[, keep_cols, drop = FALSE])
    rownames(subject_meta) <- NULL
  }

  ## ---- 5. Optionally extract gene expression --------------------------------
  gene_expr <- NULL
  if (!is.null(gene_list) && !is.null(seuratObj)) {
    if (requireNamespace("Seurat", quietly = TRUE)) {
      expr_mat <- Seurat::GetAssayData(seuratObj, slot = "data")
      genes_present <- intersect(gene_list, rownames(expr_mat))
      if (length(genes_present) > 0) {
        # Transpose to cells x genes, attach unit label, aggregate by unit
        ge_df <- as.data.frame(t(as.matrix(expr_mat[genes_present, , drop = FALSE])))
        ge_df[[unit_var]] <- meta[[unit_var]]
        agg_formula <- stats::as.formula(paste(". ~", unit_var))
        gene_expr <- stats::aggregate(
          agg_formula,
          data = ge_df,
          FUN  = mean
        )
      } else {
        warning("None of the requested genes were found in the Seurat assay.")
      }
    }
  }

  ## ---- 6. Assemble result --------------------------------------------------
  result <- list(
    proportions      = proportions_df,
    subject_metadata = subject_meta,
    gene_expression  = gene_expr,
    params           = list(
      cluster_var = cluster_var,
      unit_var    = unit_var,
      subject_var = subject_var,
      disease_var = disease_var
    )
  )
  class(result) <- "structureR_data"
  result
}

#' @export
print.structureR_data <- function(x, ...) {
  cat("structureR Data Object\n")
  cat("----------------------\n")
  cat("Units (", x$params$unit_var, "):", nrow(x$proportions), "\n")
  cat("Cluster proportions: ", ncol(x$proportions) - 1L, "clusters\n")
  if (!is.null(x$subject_metadata)) {
    cat("Subject metadata columns:", ncol(x$subject_metadata), "\n")
  }
  if (!is.null(x$gene_expression)) {
    n_genes <- ncol(x$gene_expression) - 1L
    cat("Gene expression variables:", n_genes, "\n")
  }
  invisible(x)
}
