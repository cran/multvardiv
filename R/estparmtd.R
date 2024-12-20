estparmtd <- function(x, eps = 1e-6, display = FALSE, plot = display) {
  #' Estimation of the Parameters of a Multivariate \eqn{t} Distribution
  #'
  #' Estimation of the degrees of freedom, mean vector and correlation matrix of a multivariate \eqn{t} distribution (MTD).
  #'
  #' @aliases estparmtd
  #'
  #' @usage estparmtd(x, eps = 1e-6, display = FALSE, plot = display)
  #' @param x numeric matrix or data frame.
  #' @param eps numeric. Precision for the estimation of the parameters.
  #' @param display logical. When \code{TRUE} the value of the \code{nu} parameter at each iteration is printed.
  #' @param plot logical. When \code{TRUE} the successive values of the \code{nu} parameter are plotted, allowing to visualise its convergence.
  #' @return A list of 3 elements:
  #' \itemize{
  #' \item \code{nu} non-negative numeric value. The degrees of freedom.
  #' \item \code{mu} the mean vector.
  #' \item \code{Sigma}: symmetric positive-definite matrix. The correlation matrix.
  #' }
  #' with two attributes \code{attr(, "epsilon")} (precision of the result) and \code{attr(, "k")} (number of iterations).
  #' 
  #' @details The EM method is used to estimate the parameters.
  #'
  #' @author Pierre Santagostini, Nizar Bouhlel
  #' 
  #' @references Doğru, F., Bulut, Y. M. and Arslan, O. (2018).
  #' Doubly reweighted estimators for the parameters of the multivariate t-distribution.
  #' Communications in Statistics - Theory and Methods. 47.
  #' \doi{10.1080/03610926.2018.1445861}.
  #'
  #' @seealso \code{\link{dmtd}}: probability density of a MTD
  #' 
  #' \code{\link{rmtd}}: random generation from a MTD.
  #' 
  #' @examples
  #' nu <- 3
  #' mu <- c(0, 1, 4)
  #' Sigma <- matrix(c(1, 0.6, 0.2, 0.6, 1, 0.3, 0.2, 0.3, 1), nrow = 3)
  #' x <- rmtd(100, nu, mu, Sigma)
  #' 
  #' # Estimation of the parameters
  #' estparmggd(x)
  #'
  #' @importFrom graphics plot
  #' @export
  
  if (is.vector(x))
    x <- cbind(x)
  if (is.data.frame(x))
    x <- as.matrix(x)
  
  # Number of variables
  p <- ncol(x)
  # Number of observations
  n <- nrow(x)
  
  # Parameter initialization
  nu <- 1
  mu <- rep(0, p)
  Sigma <- diag(1, nrow = p, ncol = p)
  
  k <- 0
  nu1 <- Inf
  if (plot) nuseq <- nu
  if (display) cat(nu, "  ")
  
  while (abs(nu - nu1) > eps) {
    k <- k + 1
    
    x0 <- x - matrix(mu, nrow = n, ncol = p, byrow = TRUE)
    
    invSigma <- x0 %*% solve(Sigma)
    y <- rowSums(invSigma * x0)
    
    E1 <- (nu + p) / (nu + y)
    matE1 <- matrix(E1, nrow = n, ncol = p, byrow = FALSE)
    E2 <- digamma((nu + p)/2) - log(0.5*(nu + y))
    
    # Update mu
    mu <- colSums(x * matE1) / sum(E1)
    
    # Update Sigma
    Sigma <- matrix(0, nrow = p, ncol = p)
    for (i in 1:n) {
      Sigma <- Sigma + E1[i] * cbind(x0[i, ]) %*% rbind(x0[i, ])
    }
    Sigma <- Sigma/n
    
    # Update nu
    nu1 <- nu
    nu <- .estnu(E1, E2, nu, p)
    
    if (plot) nuseq <- c(nuseq, nu)
    if (display) cat(nu, "  ")
  }
  
  if (plot) plot(nuseq, type = "l")
  
  # Returned result
  result <- list(nu = nu, mu = mu, Sigma = Sigma)
  
  # Attributes:
  # precision of the result
  attr(result, "epsilon") <- eps
  # number of iterations
  attr(result, "k") <- k
  
  return(result)
}

.estnu <- function(E1, E2, nu0, p, eps = .Machine$double.eps) {
  #' @importFrom stats uniroot
  
  # Search the zero of the function in ]0; 2*nu0]
  N <- length(E1)
  fnu <- function(z) {
    1 + log(z/2) - digamma(z/2) - sum(E1 - E2)/N
  }
  result <- uniroot(fnu, c(eps, 2*ceiling(nu0)))$root
  return(result)
}
