# Dimension p = 1

Sigma1 <- 0.5
Sigma2 <- 1

kl1_12 <- kld(Sigma1, Sigma2, distribution = "mtd", nu1 = 1, nu2 = 1, eps = 1e-16)
kl1_21 <- kld(Sigma2, Sigma1, distribution = "mtd", nu1 = 1, nu2 = 1, eps = 1e-16)

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

# Dimension p = 1, 2nd example

Sigma1 <- 0.5
nu1 <- 2
Sigma2 <- 1
nu2 <- 1

kl1_2 <- kld(Sigma1, Sigma2, distribution = "mtd", nu1 = nu1, nu2 = nu2, eps = .Machine$double.eps)

test_that("kl1_2 works (dim 1, 2nd exple)", {
  expect_equal(
    attr(kl1_2, "epsilon"), .Machine$double.eps
  )
  expect_equal(
    round(as.numeric(kl1_2), 16), 0.1447298858494002
  )
})

# Dimension p = 1, lambda*nu1/nu2 == 1
nu1 <- 2; Sigma1 <- 0.5
nu2 <- 2; Sigma2 <- 0.5

kl1_12_0 <- kld(Sigma1, Sigma2, distribution = "mtd", nu1 = nu1, nu2 = nu2, eps = 1e-16)
kl1_21_0 <- kld(Sigma2, Sigma1, distribution = "mtd", nu1 = nu1, nu2 = nu2, eps = 1e-16)

test_that("kl works (dim 1, lambda*nu1/nu2 == 1)", {
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

kl2_12 <- kld(Sigma1, Sigma2, distribution = "mtd", nu1 = 1, nu2 = 1, eps = 1e-16)
kl2_21 <- kld(Sigma2, Sigma1, distribution = "mtd", nu1 = 1, nu2 = 1, eps = 1e-16)

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

kl2 <- kld(Sigma1, Sigma2, distribution = "mtd", nu1 = 1, nu2 = 1, eps = 1e-16)

test_that("kl works (dim 2, one of the eigenvalues = 1)", {
  expect_equal(
    round(as.numeric(kl2), 15),
    log(lambda) - 3/2 * 1/sqrt(1-lambda) * log((1 - sqrt(1-lambda))/(1 + sqrt(1-lambda))) - 3
  )
})

# Dimension p = 2: third example

nu1 <- 2
Sigma1 <- diag(0.5, nrow = 2)
nu2 <- 1
Sigma2 <- diag(1, nrow = 2)

lambda <- 0.5

kl2_3 <- kld(Sigma1, Sigma2, distribution = "mtd", nu1 = nu1, nu2 = nu2, eps = .Machine$double.eps)

test_that("kl2_3 works (dim 2, 2nd example)", {
  expect_equal(
    attr(kl2_3, "epsilon"), .Machine$double.eps
  )
  expect_equal(
    round(as.numeric(kl2_3), 16), 0.1931471805599454
  )
})

# Dimension p = 2, lambda*nu1/nu2 == 1
nu1 <- 2; Sigma1 <- diag(0.5, 2)
nu2 <- 2; Sigma2 <- diag(0.5, 2)

kl2_12_0 <- kld(Sigma1, Sigma2, distribution = "mtd", nu1 = nu1, nu2 = nu2, eps = 1e-16)
kl2_21_0 <- kld(Sigma2, Sigma1, distribution = "mtd", nu1 = nu1, nu2 = nu2, eps = 1e-16)

test_that("kl works (dim 2, lambda*nu1/nu2 == 1)", {
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
nu1 <- 2; nu2 <- 4
Sigma1 <- 4*rbind(c(1, 0.6, 0.2), c(0.6, 1, 0.3), c(0.2, 0.3, 1))
Sigma2 <- rbind(c(1, 0.3, 0.1), c(0.3, 1, 0.4), c(0.1, 0.4, 1))

lambda <- 0.5

kl3_12 <- kld(Sigma1, Sigma2, distribution = "mtd", nu1 = nu1, nu2 = nu2, eps = 1e-8)
kl3_21 <- kld(Sigma2, Sigma1, distribution = "mtd", nu1 = nu1, nu2 = nu2, eps = 5e-5)
test_that("kl works (dim 3)", {
  expect_equal(
    attr(kl3_12, "epsilon"), 1e-8
  )
  expect_equal(
    round(as.numeric(kl3_12), 16), 0.9297752865860369
  )
  
  expect_equal(
    attr(kl3_21, "epsilon"), 5e-5
  )
  expect_equal(
    round(as.numeric(kl3_21), 16), 0.4074954441658625
  )
})

# Dimension p = 3, 2nd example
nu1 <- 2; nu2 <- 4
Sigma1 <- 2*rbind(c(1, 0.6, 0.2), c(0.6, 1, 0.3), c(0.2, 0.3, 1))
Sigma2 <- rbind(c(1, 0.3, 0.1), c(0.3, 1, 0.4), c(0.1, 0.4, 1))

kl3 <- kld(Sigma1, Sigma2, distribution = "mtd", nu1 = nu1, nu2 = nu2, eps = 1e-16)

test_that("kl12 works (dim 3, 2nd)", {
  expect_equal(
    attr(kl3, "epsilon"), 1e-16
  )
  expect_equal(
    round(as.numeric(kl3), 16), 0.3979439491689158
  )
})

# Dimension p = 3, 3rd example
nu1 <- 2; nu2 <- 1
Sigma1 <- diag(0.5, nrow = 3)
Sigma2 <- Sigma2 <- diag(nrow = 3)

kl3_3 <- kld(Sigma1, Sigma2, distribution = "mtd", nu1 = nu1, nu2 = nu2, eps = .Machine$double.eps)

test_that("kl3_3 works (dim 3, 3rd)", {
  expect_equal(
    attr(kl3_3, "epsilon"), .Machine$double.eps
  )
  expect_equal(
    round(as.numeric(kl3_3), 16), 0.2168616606242311
  )
})

# Dimension p = 3, lambda*nu1/nu2 == 1
nu1 <- 2; Sigma1 <- diag(0.5, 3)
nu2 <- 2; Sigma2 <- diag(0.5, 3)

kl3_12_0 <- kld(Sigma1, Sigma2, distribution = "mtd", nu1 = nu1, nu2 = nu2, eps = 1e-16)
kl3_21_0 <- kld(Sigma2, Sigma1, distribution = "mtd", nu1 = nu1, nu2 = nu2, eps = 1e-16)

test_that("kl works (dim 3, lambda*nu1/nu2 == 1)", {
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
nu1 <- 2; nu2 <- 4
Sigma1 <- 4*rbind(c(1, 0.6, 0.2, 0),
                  c(0.6, 1, 0.3, 0),
                  c(0.2, 0.3, 1, 0),
                  c(0, 0, 0, 1))
Sigma2 <- rbind(c(1, 0.3, 0.1, 0),
                c(0.3, 1, 0.4, 0),
                c(0.1, 0.4, 1, 0),
                c(0, 0, 0, 1))

kl4_12 <- kld(Sigma1, Sigma2, distribution = "mtd", nu1 = nu1, nu2 = nu2, eps = 1e-4)
kl4_21 <- kld(Sigma2, Sigma1, distribution = "mtd", nu1 = nu1, nu2 = nu2, eps = 1e-4)

test_that("kl12 works (dim 4)", {
  expect_equal(
    attr(kl4_12, "epsilon"), 1e-04
  )
  expect_equal(
    round(as.numeric(kl4_12), 16), 1.039921364830607
  )

  expect_equal(
    attr(kl4_21, "epsilon"), 1e-04
  )
  expect_equal(
    round(as.numeric(kl4_21), 16), 0.5360139882566615
  )
})

# Dimension p = 4, lambda*nu1/nu2 == 1
nu1 <- 2; Sigma1 <- diag(0.5, 4)
nu2 <- 2; Sigma2 <- diag(0.5, 4)

kl4_12_0 <- kld(Sigma1, Sigma2, distribution = "mtd", nu1 = nu1, nu2 = nu2, eps = 1e-16)
kl4_21_0 <- kld(Sigma2, Sigma1, distribution = "mtd", nu1 = nu1, nu2 = nu2, eps = 1e-16)

test_that("kl works (dim 3, lambda*nu1/nu2 == 1)", {
  expect_equal(
    as.numeric(kl4_12_0),
    0
  )
  expect_equal(
    as.numeric(kl4_21_0),
    0
  )
})

# Dimension p = 3: particular case

nu1 <- 1
Sigma1 <- diag(0.5, nrow = 3)
nu2 <- 2
Sigma2 <- diag(1, nrow = 3)

lambda <- 0.5

kl3pc <- kld(Sigma1, Sigma2, distribution = "mtd", nu1 = nu1, nu2 = nu2, eps = .Machine$double.eps)

d <- -2*log((nu2 + sqrt(nu2*lambda))/(2*nu2)) + (sqrt(nu2) - sqrt(lambda))/(sqrt(nu2) + sqrt(lambda))
d <- log( (gamma((nu1+3)/2) * gamma(nu2/2) * nu2^1.5) / (gamma((nu2+3)/2) * gamma(nu1/2) * nu1^1.5) ) +
  (nu2 - nu1)/2 * ( digamma((nu1+3)/2) - digamma(nu1/2) ) - 0.5*3*(log(lambda)) -
  (nu2 + 3)/2 * d

test_that("kl3pc works (dim 3, particular case)", {
  expect_equal(
    attr(kl3pc, "epsilon"), .Machine$double.eps
  )
  expect_equal(
    round(as.numeric(kl3pc), 16), d
  )
})

# Dimension p = 3: particular case

nu1 <- 1
Sigma1 <- diag(0.5, nrow = 3)
nu2 <- 2
Sigma2 <- diag(1, nrow = 3)

lambda <- 0.5

kl3pc <- kld(Sigma1, Sigma2, distribution = "mtd", nu1 = nu1, nu2 = nu2, eps = .Machine$double.eps)

d <- -2*log((nu2 + sqrt(nu2*lambda))/(2*nu2)) + (sqrt(nu2) - sqrt(lambda))/(sqrt(nu2) + sqrt(lambda))
d <- log( (gamma((nu1+3)/2) * gamma(nu2/2) * nu2^1.5) / (gamma((nu2+3)/2) * gamma(nu1/2) * nu1^1.5) ) +
  (nu2 - nu1)/2 * ( digamma((nu1+3)/2) - digamma(nu1/2) ) - 0.5*3*(log(lambda)) -
  (nu2 + 3)/2 * d

test_that("kl3pc works (dim 3, particular case)", {
  expect_equal(
    attr(kl3pc, "epsilon"), .Machine$double.eps
  )
  expect_equal(
    round(as.numeric(kl3pc), 16), d
  )
})

# Dimension p = 4: particular case

p <- 4
nu1 <- 2
Sigma1 <- diag(0.5, nrow = p)
nu2 <- 4
Sigma2 <- diag(1, nrow = p)

lambda <- 0.5

kl4pc <- kld(Sigma1, Sigma2, distribution = "mtd", nu1 = nu1, nu2 = nu2, eps = 1e-6)

d <- nu2/(nu2 - 2*lambda) - log(2*lambda/nu2) + nu2^2*log(2*lambda/nu2) / ((nu2 - 2*lambda)^2) + 0.5
d <- log( (gamma((nu1+p)/2) * gamma(nu2/2) * nu2^(p/2)) / (gamma((nu2+p)/2) * gamma(nu1/2) * nu1^(p/2)) ) +
  (nu2 - nu1)/2 * ( digamma((nu1+p)/2) - digamma(nu1/2) ) - 0.5*p*(log(lambda)) -
  (nu2 + p)/2 * d

test_that("kl4pc works (dim 4, particular case)", {
  expect_equal(
    attr(kl4pc, "epsilon"), 1e-6
  )
  expect_equal(
    round(as.numeric(kl4pc), 6), round(d, 6)
  )
})
