test_that("PrepVineData errors without seuratObj or metadata", {
  expect_error(PrepVineData(), "Either 'seuratObj' or 'metadata' must be provided")
})

test_that("PrepVineData errors on missing required columns", {
  meta <- data.frame(foo = 1:5, bar = 1:5)
  expect_error(
    PrepVineData(metadata = meta,
                 cluster_var = "seurat_clusters",
                 unit_var    = "cDNA_ID"),
    "missing"
  )
})

test_that("PrepVineData returns a structureR_data object from metadata", {
  set.seed(42)
  meta <- data.frame(
    cDNA_ID         = rep(paste0("s", 1:5), each = 40),
    subject_id      = rep(paste0("subj", 1:5), each = 40),
    disease_stage   = rep(c("healthy", "healthy", "disease", "disease", "healthy"), each = 40),
    seurat_clusters = sample(0:3, 200, replace = TRUE),
    age             = rep(c(25, 30, 45, 50, 35), each = 40),
    stringsAsFactors = FALSE
  )
  vd <- PrepVineData(
    metadata    = meta,
    cluster_var = "seurat_clusters",
    unit_var    = "cDNA_ID",
    subject_var = "subject_id",
    disease_var = "disease_stage"
  )

  expect_s3_class(vd, "structureR_data")
  expect_named(vd, c("proportions", "subject_metadata", "gene_expression", "params"))
})

test_that("PrepVineData proportions sum to 1 per unit", {
  set.seed(7)
  meta <- data.frame(
    cDNA_ID         = rep(paste0("u", 1:6), each = 30),
    subject_id      = rep(paste0("p", 1:6), each = 30),
    seurat_clusters = sample(0:4, 180, replace = TRUE),
    stringsAsFactors = FALSE
  )
  vd <- PrepVineData(
    metadata    = meta,
    cluster_var = "seurat_clusters",
    unit_var    = "cDNA_ID",
    subject_var = "subject_id"
  )

  prop_cols <- grep("^cluster_", colnames(vd$proportions), value = TRUE)
  row_sums  <- rowSums(vd$proportions[, prop_cols])
  expect_true(all(abs(row_sums - 1) < 1e-10))
})

test_that("PrepVineData returns correct number of rows", {
  set.seed(3)
  n_units <- 8L
  meta <- data.frame(
    cDNA_ID         = rep(paste0("u", seq_len(n_units)), each = 25),
    seurat_clusters = sample(0:2, n_units * 25, replace = TRUE),
    stringsAsFactors = FALSE
  )
  vd <- PrepVineData(metadata = meta,
                     cluster_var = "seurat_clusters",
                     unit_var    = "cDNA_ID",
                     subject_var = "missing_col")
  expect_equal(nrow(vd$proportions), n_units)
  expect_null(vd$subject_metadata)  # subject_var not found
})

test_that("PrepVineData prints without error", {
  meta <- data.frame(
    cDNA_ID         = rep("s1", 20),
    seurat_clusters = sample(0:1, 20, replace = TRUE)
  )
  vd <- PrepVineData(metadata = meta)
  expect_output(print(vd), "structureR Data Object")
})

test_that("PrepVineData requires Seurat for seuratObj input", {
  # This test only runs when Seurat is NOT installed
  skip_if(requireNamespace("Seurat", quietly = TRUE),
          "Seurat is installed; skipping absence test")
  expect_error(
    PrepVineData(seuratObj = list()),   # fake Seurat object
    "Seurat"
  )
})
