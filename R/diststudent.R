diststudent <- function(nu1, Sigma1, nu2, Sigma2, dist = c("renyi", "bhattacharyya", "hellinger"), bet = NULL, eps = 1e-06) {
  #' Distance/Divergence between Centered Multivariate \eqn{t} Distributions
  #'
  #' Computes the distance or divergence (Renyi divergence, Bhattacharyya
  #' distance or Hellinger distance) between two random vectors distributed
  #' according to multivariate \eqn{t} distributions (MTD) with zero mean vector.
  #'
  #' @aliases diststudent
  #'
  #' @usage diststudent(nu1, Sigma1, nu2, Sigma2,
  #'                    dist = c("renyi", "bhattacharyya", "hellinger"),
  #'                    bet = NULL, eps = 1e-06)
  #' @param nu1 numeric. The degrees of freedom of the first distribution.
  #' @param Sigma1 symmetric, positive-definite matrix. The correlation matrix of the first distribution.
  #' @param nu2 numeric. The degrees of freedom of the second distribution.
  #' @param Sigma2 symmetric, positive-definite matrix. The correlation matrix of the second distribution.
  #' @param dist character. The distance or divergence used.
  #' One of \code{"renyi"} (default), \code{"battacharyya"} or \code{"hellinger"}.
  #' @param bet numeric, positive and not equal to 1. Order of the Renyi divergence.
  #' Ignored if \code{distance="bhattacharyya"} or \code{distance="hellinger"}.
  #' @param eps numeric. Precision for the computation of the partial derivative of the Lauricella \eqn{D}-hypergeometric function (see Details). Default: 1e-06.
  #' @return A numeric value: the divergence between the two distributions,
  #' with two attributes \code{attr(, "epsilon")} (precision of the result of the Lauricella \eqn{D}-hypergeometric function,see Details)
  #' and \code{attr(, "k")} (number of iterations).
  #'
  #' @details Given \eqn{X_1}, a random vector of \eqn{\mathbb{R}^p} distributed according to the MTD
  #' with parameters \eqn{(\nu_1, \mathbf{0}, \Sigma_1)}
  #' and \eqn{X_2}, a random vector of \eqn{\mathbb{R}^p} distributed according to the MTD
  #' with parameters \eqn{(\nu_2, \mathbf{0}, \Sigma_2)}.
  #'
  #' Let \eqn{\delta_1 = \frac{\nu_1 + p}{2} \beta}, \eqn{\delta_2 = \frac{\nu_2 + p}{2} (1 - \beta)}
  #' and \eqn{\lambda_1, \dots, \lambda_p} the eigenvalues of the square matrix \eqn{\Sigma_1 \Sigma_2^{-1}}
  #' sorted in increasing order: \deqn{\lambda_1 < \dots < \lambda_{p-1} < \lambda_p}
  #' The Renyi divergence between \eqn{X_1} and \eqn{X_2} is:
  #' \deqn{
  #' \begin{aligned}
  #' D_R^\beta(\mathbf{X}_1||\mathbf{X}_1) &  = & \displaystyle{\frac{1}{\beta - 1} \bigg[ \beta \ln\left(\frac{\Gamma\left(\frac{\nu_1+p}{2}\right) \Gamma\left(\frac{\nu_2}{2}\right) \nu_2^{\frac{p}{2}}}{\Gamma\left(\frac{\nu_2+p}{2}\right) \Gamma\left(\frac{\nu_1}{2}\right) \nu_1^{\frac{p}{2}}}\right) + \ln\left(\frac{\Gamma\left(\frac{\nu_2+p}{2}\right)}{\Gamma\left(\frac{\nu_2}{2}\right)}\right) + \ln\left(\frac{\Gamma\left(\delta_1 + \delta_2 - \frac{p}{2}\right)}{\Gamma(\delta_1 + \delta_2)}\right) } \\
  #' && \displaystyle{- \frac{\beta}{2} \sum_{i=1}^p{\ln\lambda_i} + \ln F_D \bigg]}
  #' \end{aligned}
  #' }
  #' with \eqn{F_D} given by:
  #' \itemize{
  #' \item If \eqn{\displaystyle{\frac{\nu_1}{\nu_2} \lambda_1 > 1}}:
  #' 
  #' \eqn{
  #' \displaystyle{ F_D = F_D^{(p)}{\bigg( \delta_1, \underbrace{\frac{1}{2}, \dots, \frac{1}{2}}_p; \delta_1+\delta_2; 1-\frac{\nu_2}{\nu_1 \lambda_1}, \dots, 1-\frac{\nu_2}{\nu_1 \lambda_p} \bigg)} }
  #' }
  #' \item If \eqn{\displaystyle{\frac{\nu_1}{\nu_2} \lambda_p < 1}}:
  #' 
  #' \eqn{
  #' \displaystyle{ F_D = \prod_{i=1}^p{\left(\frac{\nu_1}{\nu_2} \lambda_i\right)^{\frac{1}{2}}} F_D^{(p)}\bigg(\delta_2, \underbrace{\frac{1}{2}, \dots, \frac{1}{2}}_p; \delta_1+\delta_2; 1-\frac{\nu_1}{\nu_2}\lambda_1, \dots, 1-\frac{\nu_1}{\nu_2}\lambda_p\bigg) }
  #' }
  #' \item If \eqn{\displaystyle{\frac{\nu_1}{\nu_2} \lambda_1 < 1}} and \eqn{\displaystyle{\frac{\nu_1}{\nu_2} \lambda_p > 1}}:
  #'
  #' \eqn{
  #' \displaystyle{ F_D = \left(\frac{\nu_2}{\nu_1} \frac{1}{\lambda_p}\right)^{\delta_2} \prod_{i=1}^p\left(\frac{\nu_1}{\nu_2}\lambda_i\right)^\frac{1}{2} F_D^{(p)}\bigg(\delta_2, \underbrace{\frac{1}{2}, \dots, \frac{1}{2}}_p, \delta_1+\delta_2-\frac{p}{2}; \delta_1+\delta2; 1-\frac{\lambda_1}{\lambda_p}, \dots, 1-\frac{\lambda_{p-1}}{\lambda_p}, 1-\frac{\nu_2}{\nu_1}\frac{1}{\lambda_p}\bigg) }
  #' }
  #' }
  #'
  #' where \eqn{F_D^{(p)}} is the Lauricella \eqn{D}-hypergeometric function defined for \eqn{p} variables:
  #' \deqn{ \displaystyle{ F_D^{(p)}\left(a; b_1, ..., b_p; g; x_1, ..., x_p\right) = \sum\limits_{m_1 \geq 0} ... \sum\limits_{m_p \geq 0}{ \frac{ (a)_{m_1+...+m_p}(b_1)_{m_1} ... (b_p)_{m_p} }{ (g)_{m_1+...+m_p} } \frac{x_1^{m_1}}{m_1!} ... \frac{x_p^{m_p}}{m_p!} } } }
  #' Its computation uses the \code{\link{lauricella}} function.
  #'
  #' The Bhattacharyya distance is given by:
  #' \deqn{D_B(\mathbf{X}_1||\mathbf{X}_2) = \frac{1}{2} D_R^{1/2}(\mathbf{X}_1||\mathbf{X}_2)}
  #'
  #' And the Hellinger distance is given by:
  #' \deqn{D_H(\mathbf{X}_1||\mathbf{X}_2) = 1 - \exp{\left(-\frac{1}{2} D_R^{1/2}(\mathbf{X}_1||\mathbf{X}_2)\right)}}
  #'
  #' @author Pierre Santagostini, Nizar Bouhlel
  #' @references N. Bouhlel and D. Rousseau (2023), Exact Rényi and Kullback-Leibler Divergences Between Multivariate t-Distributions, IEEE Signal Processing Letters.
  #' \doi{10.1109/LSP.2023.3324594}
  #'
  #' @examples
  #' nu1 <- 2
  #' Sigma1 <- matrix(c(2, 1.2, 0.4, 1.2, 2, 0.6, 0.4, 0.6, 2), nrow = 3)
  #' nu2 <- 4
  #' Sigma2 <- matrix(c(1, 0.3, 0.1, 0.3, 1, 0.4, 0.1, 0.4, 1), nrow = 3)
  #'
  #' # Renyi divergence
  #' diststudent(nu1, Sigma1, nu2, Sigma2, bet = 0.25)
  #' diststudent(nu2, Sigma2, nu1, Sigma1, bet = 0.25)
  #'
  #' # Bhattacharyya distance
  #' diststudent(nu1, Sigma1, nu2, Sigma2, dist = "bhattacharyya")
  #' diststudent(nu2, Sigma2, nu1, Sigma1, dist = "bhattacharyya")
  #'
  #' # Hellinger distance
  #' diststudent(nu1, Sigma1, nu2, Sigma2, dist = "hellinger")
  #' diststudent(nu2, Sigma2, nu1, Sigma1, dist = "hellinger")
  #'
  #' @export

  # Distances/divergences available: "renyi", "bhattacharyya" or "hellinger"
  dist <- match.arg(dist)

  switch(dist,
         renyi = {
           if (is.null(bet))
             stop("The bet argument must be provided for the Renyi distance.")

           # We must have: bet>1 and bet!=1
           if (bet <= 0 | bet == 1)
             stop("bet must be positive, different to 1.")
         },
         bhattacharyya = {
           if (!is.null(bet))
             message("bet is omitted when `dist='bhattacharyya'.")
           bet <- 0.5
         },
         hellinger = {
           if (!is.null(bet))
             message("bet is omitted when `dist='bhattacharyya'.")
           bet <- 0.5
         }
  )

  # Sigma1 and Sigma2 must be matrices
  if (is.numeric(Sigma1) & !is.matrix(Sigma1))
    Sigma1 <- matrix(Sigma1)
  if (is.numeric(Sigma2) & !is.matrix(Sigma2))
    Sigma2 <- matrix(Sigma2)

  # Number of variables
  p <- nrow(Sigma1)

  # Sigma1 and Sigma2 must be square matrices with the same size
  if (ncol(Sigma1) != p | nrow(Sigma2) != p | ncol(Sigma2) != p)
    stop("Sigma1 et Sigma2 must be square matrices with rank p.")

  # IS Sigma1 symmetric, positive-definite?
  if (!isSymmetric(Sigma1))
    stop("Sigma1 must be a symmetric, positive-definite matrix.")
  lambda1 <- eigen(Sigma1, only.values = TRUE)$values
  if (any(lambda1 < .Machine$double.eps))
    stop("Sigma1 must be a symmetric, positive-definite matrix.")

  # Is Sigma2 symmetric, positive-definite?
  if (!isSymmetric(Sigma2))
    stop("Sigma2 must be a symmetric, positive-definite matrix.")
  lambda2 <- eigen(Sigma2, only.values = TRUE)$values
  if (any(lambda2 < .Machine$double.eps))
    stop("Sigma2 must be a symmetric, positive-definite matrix.")

  # Eigenvalues of Sigma1 %*% inv(Sigma2)
  lambda <- sort(eigen(Sigma1 %*% solve(Sigma2), only.values = TRUE)$values, decreasing = FALSE)
  loglambda <- log(lambda)
  nulambda <- lambda*nu1/nu2

  delta1 <- bet*(nu1 + p)/2; delta2 <- (1 - bet)*(nu2 + p)/2

  if (nulambda[1] > 1) {
    lauric <- lauricella(delta1, rep(0.5, p), delta1 + delta2, 1 - 1/nulambda, eps = eps)
  } else if (nulambda[p] < 1) {
    lauric <- lauricella(delta2, rep(0.5, p), delta1 + delta2, 1 - nulambda, eps = eps)
    lauric <- prod(sqrt(nulambda)) * lauric
  } else {
    lauric <- lauricella(delta2, c(rep(0.5, p-1), delta1 + delta2 - p/2),
                         delta1 + delta2,
                         c(1 - lambda[1:(p-1)]/lambda[p], 1 - 1/nulambda[p]),
                         eps = eps)
    lauric <- (nulambda[p])^(-delta2) * prod(sqrt(nulambda)) * lauric
  }

  gamma1p <- gamma((nu1 + p)/2); gamma2p <- gamma((nu2 + p)/2)
  gamma1 <- gamma(nu1/2); gamma2 <- gamma(nu2/2)

  renyi <- (
    bet*log((gamma1p*gamma2*nu2^(p/2))/(gamma2p*gamma1*nu1^(p/2))) +
      log(gamma2p/gamma2) + log(gamma(delta1 + delta2 - p/2)/gamma(delta1 + delta2)) -
      bet/2 * sum(loglambda) + log(lauric)
  )/(bet - 1)

  switch(dist,
         renyi =
           return(renyi),
         bhattacharyya =
           return(renyi/2),
         hellinger =
           return(1 - exp(-renyi/2))
  )
}
