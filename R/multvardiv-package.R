#' Tools for Multivariate Probability Distributions
#'
#' This package provides tools for multivariate probability distributions:
#' \itemize{
#' \item Multivariate generalized Gaussian distribution (MGGD)
#' \item Multivariate Cauchy distribution (MCD)
#' \item Multivariate \eqn{t} distribution (MTD)
#' }
#' 
#' For each of these probability distributions, are provided:
#' \itemize{
#' \item Calculation of distances/divergences between two distributions:
#' \itemize{
#' \item Renyi divergence, Bhattacharyya distance, Hellinger distance (for now, only available for the \eqn{t} distribution): \code{\link{diststudent}}
#' \item Kullback-Leibler divergence: \code{\link{kld}}
#' }
#' \item Tools for these probability distributions: probability distribution function, parameter estimation, sample estimation
#' \itemize{
#' \item MGGD: \code{\link{dmggd}}, \code{\link{estparmggd}}, \code{\link{rmggd}}
#' \item MCD: \code{\link{dmcd}}, \code{\link{estparmcd}}, \code{\link{rmcd}}
#' \item MTD: \code{\link{dmtd}}, \code{\link{estparmtd}}, \code{\link{rmtd}}
#' \item Plot of the density of a distribution with 2 variables: \code{\link{plotmvd}}, \code{\link{contourmvd}}
#' }
#' }
#'
#' @name multvardiv-package
#' @aliases multvardiv-package multvardiv
#' @docType package
#' @author Pierre Santagostini <pierre.santagostini@institut-agro.fr>,
#' Nizar Bouhlel <nizar.bouhlel@institut-agro.fr>, <nizar.bouhlel@inrae.fr>
#' @references
#' N. Bouhlel, A. Dziri, Kullback-Leibler Divergence Between Multivariate Generalized Gaussian Distributions. IEEE Signal Processing Letters, vol. 26 no. 7, July 2019.
#' \doi{10.1109/LSP.2019.2915000}
#' 
#' N. Bouhlel, D. Rousseau, A Generic Formula and Some Special Cases for the Kullback–Leibler Divergence between Central Multivariate Cauchy Distributions. Entropy, 24, 838, July 2022.
#' \doi{10.3390/e24060838}
#'
#' N. Bouhlel and D. Rousseau (2023), Exact Rényi and Kullback-Leibler Divergences Between Multivariate t-Distributions, IEEE Signal Processing Letters.
#' \doi{10.1109/LSP.2023.3324594}
#' @keywords internal
"_PACKAGE"

NULL
