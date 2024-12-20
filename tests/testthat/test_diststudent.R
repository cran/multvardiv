nu1 <- 2
nu2 <- 4

bet <- 0.5

# Dimension p = 3

Sigma1 <- 2*rbind(c(1, 0.6, 0.2), c(0.6, 1, 0.3), c(0.2, 0.3, 1))
Sigma2 <- rbind(c(1, 0.3, 0.1), c(0.3, 1, 0.4), c(0.1, 0.4, 1))

d12_3 <- diststudent(nu1 = nu1, Sigma1 = Sigma1, nu2 = nu2, Sigma2 = Sigma2, bet = bet, eps = 1e-15)
d21_3 <- diststudent(nu1 = nu1, Sigma1 = Sigma2, nu2 = nu2, Sigma2 = Sigma1, bet = bet, eps = 1e-10)

test_that("renyi works (dim 3)", {
  expect_equal(attr(d12_3, "eps"), 1e-15)
  expect_equal(attr(d21_3, "eps"), 1e-15)
  
  expect_equal(round(as.numeric(d12_3), 15), 0.154150335789182)
  expect_equal(round(as.numeric(d21_3), 10), 0.077595179539573)
})

# Dimension p = 4

Sigma1 <- 2*rbind(c(1, 0.6, 0.2, 0), c(0.6, 1, 0.3, 0),
                  c(0.2, 0.3, 1, 0), c(0, 0, 0, 1))
Sigma2 <- rbind(c(1, 0.3, 0.1, 0), c(0.3, 1, 0.4, 0),
                c(0.1, 0.4, 1, 0), c(0, 0, 0, 1))

d12_4 <- diststudent(nu1 = nu1, Sigma1 = Sigma1, nu2 = nu2, Sigma2 = Sigma2, bet = bet, eps = 1e-6)
d21_4 <- diststudent(nu1 = nu1, Sigma1 = Sigma2, nu2 = nu2, Sigma2 = Sigma1, bet = bet, eps = 1e-4)

test_that("renyi works (dim 4)", {
  expect_equal(attr(d12_4, "eps"), 1e-6)
  expect_equal(attr(d21_4, "eps"), 1e-4)
  
  expect_equal(as.numeric(d12_4), 0.1754863591353835)
  expect_equal(as.numeric(d21_4), 0.09467178986499514)
})

# Dimension p = 4, 2nd example

Sigma1 <- 2*rbind(c(1, 0.6, 0.2, 0), c(0.6, 1, 0.3, 0),
                  c(0.2, 0.3, 1, 0), c(0, 0, 0, 1))
Sigma2 <- rbind(c(1, 0.3, 0.1, 0), c(0.3, 1, 0.4, 0),
                c(0.1, 0.4, 1, 0), c(0, 0, 0, 4))

d_4 <- diststudent(nu1 = nu1, Sigma1 = Sigma1, nu2 = nu2, Sigma2 = Sigma2, bet = bet, eps = 5e-6)

test_that("renyi works (dim 4, nd)", {
  expect_equal(attr(d_4, "eps"), 5e-6)

  expect_equal(as.numeric(d_4), 0.2223368804066479)
})
