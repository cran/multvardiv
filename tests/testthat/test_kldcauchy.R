# Dimension p = 1
Sigma1 <- 0.5
Sigma2 <- 1

kl1_12 <- kld(Sigma1, Sigma2, distribution = "mcd", eps = 1e-16)
kl1_21 <- kld(Sigma2, Sigma1, distribution = "mcd", eps = 1e-16)

lambda <- 0.5

test_that("kl works (dim 1)", {
  expect_equal(
    round(as.numeric(kl1_21), 15),
    log( (1 + sqrt(lambda))^2 / (4*sqrt(lambda)) )
  )
  expect_equal(
    round(as.numeric(kl1_21), 15),
    log( (1 + sqrt(lambda))^2 / (4*sqrt(lambda)) )
  )
})

# Dimension p = 1, Sigma1 = Sigma2
Sigma1 <- Sigma2 <- 1

kl1_12_0 <- kld(Sigma1, Sigma2, distribution = "mcd", eps = 1e-16)
kl1_21_0 <- kld(Sigma2, Sigma1, distribution = "mcd", eps = 1e-16)

test_that("kl works (dim 1, Sigma1==Sigma2)", {
  expect_equal(
    as.numeric(kl1_12_0),
    0
  )
  expect_equal(
    as.numeric(kl1_21_0),
    0
  )
})

# Dimension p = 2

Sigma1 <- diag(0.5, nrow = 2)
Sigma2 <- diag(1, nrow = 2)

kl2_12 <- kld(Sigma1, Sigma2, distribution = "mcd", eps = 1e-16)
kl2_21 <- kld(Sigma2, Sigma1, distribution = "mcd", eps = 1e-16)

lambda <- as.complex(0.5)

test_that("kl works (dim 2)", {
  expect_equal(
    round(as.numeric(kl2_12), 15),
    Re(-log(lambda) + 3/sqrt(1-1/lambda) * log(sqrt(lambda) + sqrt(lambda-1)) - 3)
  )
  expect_equal(
    round(as.numeric(kl2_21), 15),
    Re(log(lambda) + 3/sqrt(1-lambda) * log(sqrt(1/lambda) + sqrt(1/lambda-1)) - 3)
  )
})

# Dimension p = 2; 2nd example
Sigma1 <- matrix(c(0.5, 0, 0, 1), nrow = 2)
Sigma2 <- diag(nrow = 2)

lambda <- 0.5

kl2 <- kld(Sigma1, Sigma2, distribution = "mcd", eps = 1e-16)

test_that("kl works (dim 2, one of the eigenvalues = 1)", {
  expect_equal(
    round(as.numeric(kl2), 15),
    log(lambda) - 3/2 * 1/sqrt(1-lambda) * log((1 - sqrt(1-lambda))/(1 + sqrt(1-lambda))) - 3
  )
})

# Dimension p = 2, Sigma1 = Sigma2
Sigma1 <- Sigma2 <- diag(1, nrow = 2)

kl2_12_0 <- kld(Sigma1, Sigma2, distribution = "mcd", eps = 1e-16)
kl2_21_0 <- kld(Sigma2, Sigma1, distribution = "mcd", eps = 1e-16)

test_that("kl works (dim 2, Sigma1==Sigma2)", {
  expect_equal(
    as.numeric(kl2_12_0),
    0
  )
  expect_equal(
    as.numeric(kl2_21_0),
    0
  )
})

#Dimension p = 3
Sigma1 <- diag(0.5, nrow = 3)
Sigma2 <- diag(nrow = 3)

lambda <- 0.5

kl3 <- kld(Sigma1, Sigma2, distribution = "mcd", eps = 1e-16)

test_that("kl works (dim 3)", {
  expect_equal(
    round(as.numeric(kl3), 15),
    -3/2*log(lambda) + 4*log(0.5 + sqrt(lambda)/2) - 2*((1 - sqrt(lambda))/(1 + sqrt(lambda)))
  )
})

# Dimension p = 3, 2nd example
Sigma1 <- diag(c(0.5, 1, 1), nrow = 3)
Sigma2 <- diag(nrow = 3)

lambda <- 0.5

kl3 <- kld(Sigma1, Sigma2, distribution = "mcd", eps = 1e-16)

test_that("kl works (dim 3, 2 larger eigenvalues = 1)", {
  expect_equal(
    round(as.numeric(kl3), 15),
    -0.5*log(lambda) + 4*log((1 + sqrt(lambda))/2) + 2*(1 - sqrt(lambda))/(1 + sqrt(lambda))
  )
})

# Dimension p = 3, Sigma1 = Sigma2
Sigma1 <- Sigma2 <- diag(1, nrow = 3)

kl3_12_0 <- kld(Sigma1, Sigma2, distribution = "mcd", eps = 1e-16)
kl3_21_0 <- kld(Sigma2, Sigma1, distribution = "mcd", eps = 1e-16)

test_that("kl works (dim 3, Sigma1==Sigma2)", {
  expect_equal(
    as.numeric(kl3_12_0),
    0
  )
  expect_equal(
    as.numeric(kl3_21_0),
    0
  )
})

# Dimension p = 4
Sigma1 <- diag(1, 4)
Sigma2 <- matrix(c(0.5, 0, 0, 0, 0, 0.4, 0, 0, 0, 0, 0.3, 0, 0, 0, 0, 0.2), nrow = 4)

kl4.12 <- kld(Sigma1, Sigma2, distribution = "mcd", eps = 1e-3)
kl4.21 <- kld(Sigma2, Sigma1, distribution = "mcd", eps = 1e-3)

test_that("kl12 works (dim 4)", {
  expect_equal(
    round(as.numeric(kl4.12), 16),
    0.2450294025717286
  )
})

test_that("kl21 works (dim 4)", {
  expect_equal(
    round(as.numeric(kl4.21), 16),
    0.2632419419270904
  )
})

# Dimension p = 4, Sigma1 = Sigma2
Sigma1 <- Sigma2 <- diag(1, nrow = 4)

kl4_12_0 <- kld(Sigma1, Sigma2, distribution = "mcd", eps = 1e-6)
kl4_21_0 <- kld(Sigma2, Sigma1, distribution = "mcd", eps = 1e-6)

test_that("kl works (dim 4, Sigma1==Sigma2)", {
  expect_equal(
    as.numeric(kl4_12_0),
    0
  )
  expect_equal(
    as.numeric(kl4_21_0),
    0
  )
})
