# this test checks that the package builds and links correctly

test_that("Exported R functions exist", {
  expect_true(is.function(dynConfiRIDP::fit_rtconf_models_formula))
  expect_true(is.function(dynConfiRIDP::fit_rtconf_formula))
})

test_that("Cpp functions are compiled and linked", {
  # checkk if the underlying .Call symbol is loaded
  expect_true(is.function(dynConfiRIDP:::grid_search_worker))
  expect_true(is.function(dynConfiRIDP:::nlopt_optimizer))
})