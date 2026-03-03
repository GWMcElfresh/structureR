skip_if_not_installed("rvinecopulib")

make_pseudo <- function(n = 60, p = 4, seed = 42) {
  set.seed(seed)
  m <- matrix(runif(n * p), nrow = n, ncol = p)
  colnames(m) <- c("disease_stage", "prop_0", "prop_1", "age")
  m
}

test_that("FitFactorVine returns a structureR_vine object", {
  pseudo <- make_pseudo()
  fit    <- FitFactorVine(pseudo, factor_var = "disease_stage")

  expect_s3_class(fit, "structureR_vine")
  expect_named(fit, c("vine", "marginals", "tree_plot", "params"))
})

test_that("FitFactorVine factor index is computed correctly", {
  pseudo <- make_pseudo()
  fit    <- FitFactorVine(pseudo, factor_var = "prop_1")

  expect_equal(fit$marginals$factor_var, "prop_1")
  # After reordering, column 1 should be "prop_1"
  expect_equal(colnames(fit$marginals$ordered_pseudo)[1], "prop_1")
})

test_that("FitFactorVine falls back gracefully when factor_var not found", {
  pseudo <- make_pseudo()
  expect_warning(
    fit <- FitFactorVine(pseudo, factor_var = "nonexistent"),
    "not found"
  )
  expect_s3_class(fit, "structureR_vine")
})

test_that("FitFactorVine accepts numeric factor_var index", {
  pseudo <- make_pseudo()
  fit    <- FitFactorVine(pseudo, factor_var = 3L)

  expect_equal(fit$marginals$factor_idx, 3L)
})

test_that("FitFactorVine prints without error", {
  pseudo <- make_pseudo()
  fit    <- FitFactorVine(pseudo)
  expect_output(print(fit), "structureR Factor-Vine Model")
})

test_that("FitFactorVine errors on pseudo_obs outside (0,1)", {
  bad <- matrix(c(0.5, 1.5, 0.3, 0.7), nrow = 2, ncol = 2)
  expect_error(FitFactorVine(bad), "strictly in \\(0, 1\\)")
})

test_that("FitFactorVine summary runs without error", {
  pseudo <- make_pseudo()
  fit    <- FitFactorVine(pseudo)
  expect_output(summary(fit), "structureR Factor-Vine Model")
})

test_that("FitFactorVine var_names override works", {
  set.seed(10)
  pseudo <- matrix(runif(120), nrow = 30, ncol = 4)
  fit    <- FitFactorVine(pseudo, var_names = c("D", "A", "B", "C"),
                          factor_var = "D")
  expect_equal(fit$marginals$var_names, c("D", "A", "B", "C"))
})
