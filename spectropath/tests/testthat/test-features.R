testthat::test_that("path scalar features match path_features", {
  u <- c(-2, -1, 0, 0.8, 2)
  f <- c(0, 0.5, 1, 0.2, 0)
  path <- cbind(u, f)

  got <- path_features(path, depth = 4, normalize = FALSE)

  testthat::expect_equal(got$levy_area, levy_area(path, normalize = FALSE))
  testthat::expect_equal(got$curl_u, curl_u(path, normalize = FALSE))
  testthat::expect_equal(got$curl_f, curl_f(path, normalize = FALSE))
  testthat::expect_equal(got$skew_signature, skew_signature(path, normalize = FALSE))
  testthat::expect_equal(got$jerk_u, jerk_u(path, normalize = FALSE))
  testthat::expect_equal(got$twist_uf, twist_uf(path, normalize = FALSE))
  testthat::expect_equal(got$jerk_f, jerk_f(path, normalize = FALSE))
})

testthat::test_that("emission absorption area changes sign under ordering reversal", {
  u <- seq(-5, 5, length.out = 400)
  f1 <- exp(-0.5 * ((u - 1) / 0.7)^2) -
    0.8 * exp(-0.5 * ((u + 1) / 0.4)^2)
  f2 <- rev(f1)

  a1 <- emabs_area(cbind(u, f1), normalize = TRUE)
  a2 <- emabs_area(cbind(u, f2), normalize = TRUE)

  testthat::expect_equal(a1, -a2, tolerance = 1e-10)
})

testthat::test_that("classical features return expected columns", {
  u <- seq(-5, 5, length.out = 300)
  f <- exp(-0.5 * (u / 0.8)^2)
  got <- classical_features(cbind(u, f))

  testthat::expect_true(all(c(
    "line_flux", "equivalent_width", "peak_flux", "fwhm", "w80",
    "centroid", "line_dispersion", "skewness", "kurtosis_excess"
  ) %in% names(got)))
  testthat::expect_true(is.finite(got$w80))
  testthat::expect_true(is.finite(got$fwhm))
})
