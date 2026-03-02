skip_if_not_installed("rvinecopulib")

make_vine_fit <- function(n = 80, p = 4, seed = 1) {
  set.seed(seed)
  pseudo <- matrix(runif(n * p), nrow = n, ncol = p)
  colnames(pseudo) <- c("disease_stage", "prop_0", "prop_1", "age")
  FitFactorVine(pseudo, factor_var = "disease_stage")
}

test_that("SimulateCounterfactual returns a structureR_counterfactual object", {
  fit <- make_vine_fit()
  cf  <- SimulateCounterfactual(fit, factor_value = 0.8, n_sim = 30)

  expect_s3_class(cf, "structureR_counterfactual")
  expect_named(cf, c("samples", "factor_var", "factor_value", "n_sim"))
})

test_that("SimulateCounterfactual returns correct dimensions", {
  fit <- make_vine_fit()
  cf  <- SimulateCounterfactual(fit, factor_value = 0.2, n_sim = 50)

  expect_equal(nrow(cf$samples), 50L)
  expect_equal(ncol(cf$samples), fit$params$n_vars)
})

test_that("SimulateCounterfactual samples are in (0,1)", {
  fit <- make_vine_fit()
  cf  <- SimulateCounterfactual(fit, factor_value = 0.5, n_sim = 40)

  expect_true(all(cf$samples > 0 & cf$samples < 1))
})

test_that("SimulateCounterfactual errors on invalid fitted_vine input", {
  expect_error(
    SimulateCounterfactual(list(), factor_value = 0.5),
    "structureR_vine"
  )
})

test_that("SimulateCounterfactual errors on factor_value outside (0,1)", {
  fit <- make_vine_fit()
  expect_error(SimulateCounterfactual(fit, factor_value = 0),   "strictly in \\(0, 1\\)")
  expect_error(SimulateCounterfactual(fit, factor_value = 1),   "strictly in \\(0, 1\\)")
  expect_error(SimulateCounterfactual(fit, factor_value = 1.5), "strictly in \\(0, 1\\)")
})

test_that("SimulateCounterfactual prints without error", {
  fit <- make_vine_fit()
  cf  <- SimulateCounterfactual(fit, factor_value = 0.3, n_sim = 10)
  expect_output(print(cf), "structureR Counterfactual Simulation")
})

test_that("SimulateCounterfactual factor_value is stored correctly", {
  fit <- make_vine_fit()
  cf  <- SimulateCounterfactual(fit, factor_value = 0.75, n_sim = 20)
  expect_equal(cf$factor_value, 0.75)
  expect_equal(cf$n_sim, 20L)
})
