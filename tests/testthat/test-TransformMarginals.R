test_that("TransformMarginals works on a plain numeric data.frame", {
  set.seed(1)
  df <- data.frame(x = rnorm(50), y = rbeta(50, 2, 5), z = rpois(50, 1))
  pseudo <- TransformMarginals(df)

  expect_true(is.matrix(pseudo))
  expect_equal(ncol(pseudo), 3L)
  expect_equal(nrow(pseudo), 50L)
  expect_true(all(pseudo > 0 & pseudo < 1))
})

test_that("TransformMarginals drops non-numeric columns silently", {
  set.seed(2)
  df <- data.frame(
    id  = paste0("s", 1:30),
    val = rnorm(30),
    grp = sample(letters[1:3], 30, replace = TRUE),
    stringsAsFactors = FALSE
  )
  pseudo <- TransformMarginals(df)
  expect_equal(ncol(pseudo), 1L)   # only 'val' column
  expect_true(all(pseudo > 0 & pseudo < 1))
})

test_that("TransformMarginals errors on all non-numeric input", {
  df <- data.frame(a = letters[1:5], b = LETTERS[1:5], stringsAsFactors = FALSE)
  expect_error(TransformMarginals(df), "No numeric columns")
})

test_that("TransformMarginals works on a plain matrix", {
  set.seed(3)
  m <- matrix(rnorm(100), nrow = 20, ncol = 5)
  pseudo <- TransformMarginals(m)
  expect_equal(dim(pseudo), c(20L, 5L))
  expect_true(all(pseudo > 0 & pseudo < 1))
})

test_that("TransformMarginals zero-augmented path produces (0,1) output", {
  set.seed(4)
  # ~40 % zeros
  x_zi <- c(rep(0, 20), rexp(30, rate = 2))
  df   <- data.frame(gene_a = x_zi)
  pseudo <- TransformMarginals(df, zero_augmented_genes = "gene_a")
  expect_true(all(pseudo > 0 & pseudo < 1))
})

test_that("TransformMarginals works with use_rcpp = FALSE (pure R path)", {
  set.seed(5)
  df <- data.frame(a = rnorm(40), b = runif(40))
  pseudo_r   <- TransformMarginals(df, use_rcpp = FALSE)
  expect_true(all(pseudo_r > 0 & pseudo_r < 1))
})

test_that("TransformMarginals on structureR_data returns matrix", {
  set.seed(9)
  meta <- data.frame(
    cDNA_ID         = rep(paste0("s", 1:10), each = 20),
    subject_id      = rep(paste0("p", 1:10), each = 20),
    seurat_clusters = sample(0:2, 200, replace = TRUE),
    stringsAsFactors = FALSE
  )
  vd     <- PrepVineData(metadata = meta)
  pseudo <- TransformMarginals(vd)

  expect_true(is.matrix(pseudo))
  expect_true(all(pseudo > 0 & pseudo < 1))
  # Columns = number of cluster proportion columns in proportions
  prop_cols <- grep("^cluster_", colnames(vd$proportions), value = TRUE)
  expect_equal(ncol(pseudo), length(prop_cols))
})
