dmtd <- function(x, nu, mu, Sigma, tol = 1e-6) {
  #' Density of a Multivariate \eqn{t} Distribution
  #'
  #' Density of the multivariate (\eqn{p} variables) \eqn{t} distribution (MTD)
  #' with degrees of freedom \code{nu}, mean vector \code{mu} and
  #' correlation matrix \code{Sigma}.
  #'
  #' @aliases dmtd
  #'
  #' @usage dmtd(x, nu, mu, Sigma, tol = 1e-6)
  #' @param x length \eqn{p} numeric vector.
  #' @param nu numeric. The degrees of freedom.
  #' @param mu length \eqn{p} numeric vector. The mean vector.
  #' @param Sigma symmetric, positive-definite square matrix of order \eqn{p}. The correlation matrix.
  #' @param tol tolerance (relative to largest variance) for numerical lack of positive-definiteness in Sigma.
  #' @return The value of the density.
  #'
  #' @details The density function of a multivariate \eqn{t} distribution
  #' with \eqn{p} variables is given by:
  #' \deqn{ \displaystyle{ f(\mathbf{x}|\nu, \boldsymbol{\mu}, \Sigma) = \frac{\Gamma\left( \frac{\nu+p}{2} \right) |\Sigma|^{-1/2}}{\Gamma\left( \frac{\nu}{2} \right) (\nu \pi)^{p/2}} \left( 1 + \frac{1}{\nu} (\mathbf{x}-\boldsymbol{\mu})^T \Sigma^{-1} (\mathbf{x}-\boldsymbol{\mu}) \right)^{-\frac{\nu+p}{2}} } }
  #'
  #' When \eqn{p=1} (univariate case) it is the location-scale \eqn{t} distribution, with density function:
  #' \deqn{ \displaystyle{ f(x|\nu, \mu, \sigma^2) = \frac{\Gamma\left( \frac{\nu+1}{2} \right)}{\Gamma\left( \frac{\nu}{2} \right) \sqrt{\nu \pi \sigma^2}} \left(1 + \frac{(x-\mu)^2}{\nu \sigma^2}\right)^{-\frac{\nu+1}{2}} } }
  #'
  #' @author Pierre Santagostini, Nizar Bouhlel
  #'
  #' @references S. Kotz and Saralees Nadarajah (2004), Multivariate \eqn{t} Distributions and Their Applications, Cambridge University Press.
  #'
  #' @seealso \code{\link{rmtd}}: random generation from a MTD.
  #' 
  #' \code{\link{estparmtd}}: estimation of the parameters of a MTD.
  #' 
  #' \code{\link{plotmvd}}, \code{\link{contourmvd}}: plot of the probability density of a bivariate distribution.
  #' 
  #' @examples
  #' nu <- 1
  #' mu <- c(0, 1, 4)
  #' Sigma <- matrix(c(0.8, 0.3, 0.2, 0.3, 0.2, 0.1, 0.2, 0.1, 0.2), nrow = 3)
  #' dmtd(c(0, 1, 4), nu, mu, Sigma)
  #' dmtd(c(1, 2, 3), nu, mu, Sigma)
  #'
  #' # Univariate
  #' dmtd(1, 3, 0, 1)
  #' dt(1, 3)
  #'
  #' @export

  # Number of variables
  p <- length(mu)

  # Sigma must be a matrix
  if (is.numeric(Sigma) & !is.matrix(Sigma))
    Sigma <- as.matrix(Sigma)

  # x must have the same length as mu
  if (length(x) != p)
    stop(paste("x does not have", p, "elements.\n x and mu must have the same length."))

  # Sigma1 and Sigma2 must be square matrices with p rows and p columns
  if (NROW(Sigma) != p | NCOL(Sigma) != p)
    stop("Sigma must be a square matrix with size equal to length(mu).")

  if (p == 1) {
    return(as.numeric(
      gamma((nu+1)/2) / (gamma(nu/2)*sqrt(nu*pi*Sigma)) *
        (1 + (x-mu)^2 / (nu*Sigma))^(-(nu+1)/2)
    ))
  }

  # Is Sigma symmetric?
  if (!isSymmetric(Sigma))
    stop("Sigma must be a symmetric, positive-definite matrix.")

  # Eigenvalues and eigenvectors of Sigma
  eig <- eigen(Sigma, symmetric = TRUE)
  lambda <- eig$values

  # Is Sigma positive-definite?
  if (any(lambda < tol * max(abs(lambda))))
    stop("Sigma must be a symmetric, positive-definite matrix.")

  # Inverse of matrix Sigma
  invSigma <- solve(Sigma)

  xcent <- cbind(x - mu)

  # Computation of the density
  result <- gamma((nu+p)/2) / ( gamma(nu/2)*(nu*pi)^(p/2) )
  result <- result * (1 + t(xcent) %*% invSigma %*% xcent)^(-(nu+p)/2)

  return(as.numeric(result))
}
