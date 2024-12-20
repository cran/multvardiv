rmtd <- function(n, nu, mu, Sigma, tol = 1e-6) {
  #' Simulate from a Multivariate \eqn{t} Distribution
  #'
  #' Produces one or more samples from the multivariate (\eqn{p} variables) \eqn{t} distribution (MTD)
  #' with degrees of freedom \code{nu}, mean vector \code{mu} and
  #' correlation matrix \code{Sigma}.
  #'
  #' @aliases rmtd
  #'
  #' @usage rmtd(n, nu, mu, Sigma, tol = 1e-6)
  #' @param n integer. Number of observations.
  #' @param nu numeric. The degrees of freedom.
  #' @param mu length \eqn{p} numeric vector. The mean vector
  #' @param Sigma symmetric, positive-definite square matrix of order \eqn{p}. The correlation matrix.
  #' @param tol tolerance for numerical lack of positive-definiteness in Sigma (for \code{mvrnorm}, see Details).
  #' @return A matrix with \eqn{p} columns and \eqn{n} rows.
  #'
  #' @details A sample from a MTD with parameters \eqn{\nu}, \eqn{\boldsymbol{\mu}} and \eqn{\Sigma}
  #' can be generated using:
  #' \deqn{\displaystyle{\mathbf{X} = \boldsymbol{\mu} + \mathbf{Y} \sqrt{\frac{\nu}{u}}}}
  #' where \eqn{Y} is a random vector distributed among a centered Gaussian density
  #' with covariance matrix \eqn{\Sigma} (generated using \code{\link[MASS]{mvrnorm}})
  #' and \eqn{u} is distributed among a Chi-squared distribution with \eqn{\nu} degrees of freedom.
  #'
  #' @author Pierre Santagostini, Nizar Bouhlel
  #'
  #' @references S. Kotz and Saralees Nadarajah (2004), Multivariate \eqn{t} Distributions and Their Applications, Cambridge University Press.
  #'
  #' @seealso \code{\link{dmtd}}: probability density of a MTD.
  #'
  #' \code{\link{estparmtd}}: estimation of the parameters of a MTD.
  #'
  #' @examples
  #' nu <- 3
  #' mu <- c(0, 1, 4)
  #' Sigma <- matrix(c(1, 0.6, 0.2, 0.6, 1, 0.3, 0.2, 0.3, 1), nrow = 3)
  #' x <- rmtd(10000, nu, mu, Sigma)
  #' head(x)
  #' dim(x)
  #' mu; colMeans(x)
  #' nu/(nu-2)*Sigma; var(x)
  #'
  #' @importFrom MASS mvrnorm
  #' @importFrom stats rchisq
  #' @export

  # Number of variables
  p <- length(mu)

  # Sigma must be a matrix
  if (is.numeric(Sigma) & !is.matrix(Sigma))
    Sigma <- as.matrix(Sigma)

  # Sigma must be a square matricx with p rows and p columns
  if (nrow(Sigma) != p | ncol(Sigma) != p)
    stop("Sigma must be a square matrix with size equal to length(mu).")

  # IS Sigma symmetric?
  if (!isSymmetric(Sigma))
    stop("Sigma must be a symmetric, positive-definite matrix.")

  # Eigenvalues and eigenvectors of Sigma
  eig <- eigen(Sigma, symmetric = TRUE)
  lambda <- eig$values

  # Is Sigma positive-definite?
  if (any(lambda < tol * max(abs(lambda))))
    stop("Sigma must be a symmetric, positive-definite matrix.")

  # Computation of the density
  x <- rchisq(n, df = nu)
  result <- mvrnorm(n, rep(0, p), Sigma, tol = tol)
  result <- matrix(mu, nrow = n, ncol = p, byrow = TRUE) + result/sqrt(x/nu)

  return(result)
}
