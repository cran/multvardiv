#' Kullback-Leibler Divergence between Centered Multivariate Distributions
#'
#' Computes the Kullback-Leibler divergence between two random vectors distributed
#' according to centered multivariate distributions:
#' \itemize{
#' \item multivariate generalized Gaussian distribution (MGGD) with zero mean vector, using the \code{\link{kldggd}} function
#' \item multivariate Cauchy distribution (MCD) with zero location vector, using the \code{\link{kldcauchy}} function
#' \item multivariate \eqn{t} distribution (MTD) with zero mean vector, using the \code{\link{kldstudent}} function
#' }
#' One can also use one of the \code{\link{kldggd}}, \code{\link{kldcauchy}} or \code{\link{kldstudent}} functions, depending on the probability distribution.
#'
#' @aliases kld
#'
#' @usage kld(Sigma1, Sigma2, distribution = c("mggd", "mcd", "mtd"),
#'            beta1 = NULL, beta2 = NULL, nu1 = NULL, nu2 = NULL, eps = 1e-06)
#' @param Sigma1 symmetric, positive-definite matrix. The scatter matrix of the first distribution.
#' @param Sigma2 symmetric, positive-definite matrix. The scatter matrix of the second distribution.
#' @param distribution the probability distribution. It can be \code{"mggd"} (multivariate generalized Gaussian distribution) \code{"mcd"} (multivariate Cauchy) or \code{"mtd"} (multivariate \eqn{t}).
#' @param beta1,beta2 numeric. If \code{distribution = "mggd"}, the shape parameters of the first and second distributions.
#' \code{NULL} if \code{distribution} is \code{"mcd"} or \code{"mtd"}.
#' @param nu1,nu2 numeric. If \code{distribution = "mtd"}, the degrees of freedom of the first and second distributions.
#' \code{NULL} if \code{distribution} is \code{"mggd"} or \code{"mcd"}.
#' @param eps numeric.
#' Precision for the computation of the Lauricella \eqn{D}-hypergeometric function
#' if \code{distribution} is \code{"mggd"} (see \code{\link{kldggd}})
#' or of its partial derivative if \code{distribution = "mcd"} or \code{distribution = "mtd"} (see \code{\link{kldcauchy}} or \code{\link{kldstudent}}). Default: 1e-06.
#' @return A numeric value: the Kullback-Leibler divergence between the two distributions,
#' with two attributes \code{attr(, "epsilon")}
#' (precision of the Lauricella \eqn{D}-hypergeometric function or of its partial derivative)
#' and \code{attr(, "k")} (number of iterations).
#'
#' @author Pierre Santagostini, Nizar Bouhlel
#' @references N. Bouhlel, A. Dziri, Kullback-Leibler Divergence Between Multivariate Generalized Gaussian Distributions.
#' IEEE Signal Processing Letters, vol. 26 no. 7, July 2019. \doi{10.1109/LSP.2019.2915000}
#' 
#' N. Bouhlel, D. Rousseau, A Generic Formula and Some Special Cases for the Kullback–Leibler Divergence between Central Multivariate Cauchy Distributions.
#' Entropy, 24, 838, July 2022. \doi{10.3390/e24060838}
#' 
#' N. Bouhlel and D. Rousseau (2023), Exact Rényi and Kullback-Leibler Divergences Between Multivariate t-Distributions, IEEE Signal Processing Letters.
#' \doi{10.1109/LSP.2023.3324594}
#'
#' @examples
#' \donttest{
#' # Generalized Gaussian distributions
#' beta1 <- 0.74
#' beta2 <- 0.55
#' Sigma1 <- matrix(c(0.8, 0.3, 0.2, 0.3, 0.2, 0.1, 0.2, 0.1, 0.2), nrow = 3)
#' Sigma2 <- matrix(c(1, 0.3, 0.2, 0.3, 0.5, 0.1, 0.2, 0.1, 0.7), nrow = 3)
#' # Kullback-Leibler divergence
#' kl12 <- kld(Sigma1, Sigma2, "mggd", beta1 = beta1, beta2 = beta2)
#' kl21 <- kld(Sigma2, Sigma1, "mggd", beta1 = beta2, beta2 = beta1)
#' print(kl12)
#' print(kl21)
#' # Distance (symmetrized Kullback-Leibler divergence)
#' kldist <- as.numeric(kl12) + as.numeric(kl21)
#' print(kldist)
#' 
#' # Cauchy distributions
#' Sigma1 <- matrix(c(1, 0.6, 0.2, 0.6, 1, 0.3, 0.2, 0.3, 1), nrow = 3)
#' Sigma2 <- matrix(c(1, 0.3, 0.1, 0.3, 1, 0.4, 0.1, 0.4, 1), nrow = 3)
#' kld(Sigma1, Sigma2, "mcd")
#' kld(Sigma2, Sigma1, "mcd")
#' 
#' Sigma1 <- matrix(c(0.5, 0, 0, 0, 0.4, 0, 0, 0, 0.3), nrow = 3)
#' Sigma2 <- diag(1, 3)
#' # Case when all eigenvalues of Sigma1 %*% solve(Sigma2) are < 1
#' kld(Sigma1, Sigma2, "mcd")
#' # Case when all eigenvalues of Sigma1 %*% solve(Sigma2) are > 1
#' kld(Sigma2, Sigma1, "mcd")
#' 
#' # Student distributions
#' nu1 <- 2
#' Sigma1 <- matrix(c(2, 1.2, 0.4, 1.2, 2, 0.6, 0.4, 0.6, 2), nrow = 3)
#' nu2 <- 4
#' Sigma2 <- matrix(c(1, 0.3, 0.1, 0.3, 1, 0.4, 0.1, 0.4, 1), nrow = 3)
#' # Kullback-Leibler divergence
#' kld(Sigma1, Sigma2, "mtd", nu1 = nu1, nu2 = nu2)
#' kld(Sigma2, Sigma1, "mtd", nu1 = nu2, nu2 = nu1)
#' }
#'
#' @importFrom utils combn
#' @export

kld <- function(Sigma1, Sigma2, distribution = c("mggd", "mcd", "mtd"),
                beta1 = NULL, beta2 = NULL, nu1 = NULL, nu2 = NULL, eps = 1e-06) {
  
  distribution <- match.arg(distribution)
  
  return(
    switch(distribution,
           "mggd" = {
             kldggd(Sigma1 = Sigma1, beta1 = beta1, Sigma2 = Sigma2, beta2 = beta2,
                    eps = eps)
           },
           "mcd" = {
             kldcauchy(Sigma1 = Sigma1, Sigma2 = Sigma2, eps = eps)
           },
           "mtd" = {
             kldstudent(nu1 = nu1, Sigma1 = Sigma1, nu2 = nu2, Sigma2 = Sigma2, eps = eps)
           }
    )
  )
}
