set.seed(42)

nu1 <- 4
mu1 <- c(0, 0, 0)
Sigma1 <- matrix(c(0.8, 0.3, 0.2, 0.3, 0.2, 0.1, 0.2, 0.1, 0.2), nrow = 3)

x1 <- rmtd(10000, nu1, mu1, Sigma1)

set.seed(42)

nu2 <- 0.5
mu2 <- c(1, 2, 4)
Sigma2 <- matrix(c(1, 0.3, 0.2, 0.3, 0.5, 0.1, 0.2, 0.1, 0.7), nrow = 3)

x2 <- rmtd(10000, nu2, mu2, Sigma2)

# Esimation des paramètres
x1estim <- estparmtd(x1, eps = 1e-8, display = FALSE, plot = FALSE)
x2estim <- estparmtd(x2, eps = 1e-8, display = FALSE, plot = FALSE)

test_that("estparmggd works",
  {
    expect_equal(nu1, x1estim$nu, tolerance = 0.05)
    expect_equal(mu1, x1estim$mu, tolerance = 0.05)
    expect_equal(Sigma1, x1estim$Sigma, tolerance = 0.05)
    
    expect_equal(nu2, x2estim$nu, tolerance = 0.05)
    expect_equal(mu2, x2estim$mu, tolerance = 0.05)
    expect_equal(Sigma2, x2estim$Sigma, tolerance = 0.05)
  }
)
