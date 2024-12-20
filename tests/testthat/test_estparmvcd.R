set.seed(42)

mu1 <- c(0, 0, 0)
Sigma1 <- matrix(c(0.8, 0.3, 0.2, 0.3, 0.2, 0.1, 0.2, 0.1, 0.2), nrow = 3)

x1 <- rmcd(10000, mu1, Sigma1)

set.seed(42)

mu2 <- c(1, 2, 4)
Sigma2 <- matrix(c(1, 0.3, 0.2, 0.3, 0.5, 0.1, 0.2, 0.1, 0.7), nrow = 3)

x2 <- rmcd(10000, mu2, Sigma2)

# Esimation des paramètres
x1estim <- estparmcd(x1, eps = 1e-8)
x2estim <- estparmcd(x2, eps = 1e-8)

test_that("estparmcd works",
  {
    expect_equal(mu1, x1estim$mu, tolerance = 0.05)
    expect_equal(Sigma1, x1estim$Sigma, tolerance = 0.05)
    
    expect_equal(mu2, x2estim$mu, tolerance = 0.05)
    expect_equal(Sigma2, x2estim$Sigma, tolerance = 0.05)
  }
)
