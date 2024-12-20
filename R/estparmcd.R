estparmcd <- function(x, eps = 1e-6) {
  #' Estimation of the Parameters of a Multivariate Cauchy Distribution
  #'
  #' Estimation of the mean vector and correlation matrix of a multivariate Cauchy distribution (MCD).
  #'
  #' @aliases estparmcd
  #'
  #' @usage estparmcd(x, eps = 1e-6)
  #' @param x numeric matrix or data frame.
  #' @param eps numeric. Precision for the estimation of the parameters.
  #' @return A list of 2 elements:
  #' \itemize{
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
  #' @seealso \code{\link{dmcd}}: probability density of a MTD
  #' 
  #' \code{\link{rmcd}}: random generation from a MTD.
  #' 
  #' @examples
  #' mu <- c(0, 1, 4)
  #' Sigma <- matrix(c(1, 0.6, 0.2, 0.6, 1, 0.3, 0.2, 0.3, 1), nrow = 3)
  #' x <- rmcd(100, mu, Sigma)
  #' 
  #' # Estimation of the parameters
  #' estparmcd(x)
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
  mu <- rep(0, p)
  Sigma <- diag(1, nrow = p, ncol = p)
  
  k <- 0
  mu1 <- rep(Inf, p)
  # if (plot) nuseq <- nu
  # if (display) cat(nu, "  ")
  
  while (sqrt(sum((mu - mu1)^2)) > eps) {
    k <- k + 1
    
    x0 <- x - matrix(mu, nrow = n, ncol = p, byrow = TRUE)
    
    invSigma <- x0 %*% solve(Sigma)
    y <- rowSums(invSigma * x0)
    
    E1 <- (1 + p) / (1 + y)
    matE1 <- matrix(E1, nrow = n, ncol = p, byrow = FALSE)
    E2 <- digamma((1 + p)/2) - log(0.5*(1 + y))
    
    # Update mu
    mu1 <- mu
    mu <- colSums(x * matE1) / sum(E1)
    
    # Update Sigma
    Sigma <- matrix(0, nrow = p, ncol = p)
    for (i in 1:n) {
      Sigma <- Sigma + E1[i] * cbind(x0[i, ]) %*% rbind(x0[i, ])
    }
    Sigma <- Sigma/n
    
    # if (plot) nuseq <- c(nuseq, nu)
    # if (display) cat(nu, "  ")
  }
  
  # if (plot) plot(nuseq, type = "l")
  
  # Returned result
  result <- list(mu = mu, Sigma = Sigma)
  
  # Attributes:
  # precision of the result
  attr(result, "epsilon") <- eps
  # number of iterations
  attr(result, "k") <- k
  
  return(result)
}
