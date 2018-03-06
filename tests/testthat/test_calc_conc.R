
# Calculate standard concentration ----------------------------------------
context("Calculate standard concentration")

test_that("Calculation is correct",{
  expect_identical(calc_std_conc(c(0:7), 10000),
                   c(0, 10000/(4^6), 10000/(4^5), 10000/(4^4), 10000/(4^3), 10000/(4^2), 10000/(4), 10000))

  expect_identical(calc_std_conc(c(7:1, 0), 5, c(1:20)),
                   c(5/(c(1:7) %>% cumprod()), 0))

  expect_identical(calc_std_conc(c(0:7), 10000)[1], 0)

  expect_identical(calc_std_conc(c(7:1, 0), 5, c(1, 2, 2, 2, 4, 6, 6, 100000))[8], 0)
})

test_that("Numeric is returned", {
  expect_is(calc_std_conc(c(0:9), 5), "numeric")

  expect_is(calc_std_conc(c(letters[1:8]), 5), "numeric")
  expect_is(calc_std_conc(c(letters[1:7], 0), 5), "numeric")
  expect_is(calc_std_conc(c(letters[1:7], 0, 1), 5), "numeric")

  expect_is(calc_std_conc(c(7:1, 0), 5, c(1, 2, 2, 2, 4, 6, 6, 0)), "numeric")
  expect_is(calc_std_conc(c(7:1, 0), 5, c(1, 2, 2, 2, 4, 6, 6, 100000)), "numeric")
  expect_is(calc_std_conc(c(8:1), 5, c(1, 2, 2, 2, 4, 6, 6, 100000)), "numeric")
})

test_that("Errors and warnings are thrown", {
  expect_warning(calc_std_conc(c(0:6), 5))

  expect_error(calc_std_conc(c(7:1, 0), 5, c(1, 2)))
})
