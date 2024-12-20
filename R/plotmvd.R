#' Plot a Bivariate Density
#'
#' Plots the probability density of a multivariate distribution with 2 variables:
#' \itemize{
#' \item generalized Gaussian distribution (MGGD) with mean vector \code{mu}, dispersion matrix \code{Sigma} and shape parameter \code{beta}
#' \item Cauchy distribution (MCD) with location parameter \code{mu} and scatter matrix \code{Sigma}
#' \item \eqn{t} distribution (MTD) with location parameter \code{mu} and scatter matrix \code{Sigma}
#' }
#' This function uses the \code{\link[rgl]{plot3d.function}} function.
#' 
#' @aliases plotmvd
#'
#' @usage plotmvd(mu, Sigma, beta = NULL, nu = NULL,
#'                distribution = c("mggd", "mcd", "mtd"),
#'                xlim = c(mu[1] + c(-10, 10)*Sigma[1, 1]),
#'                ylim = c(mu[2] + c(-10, 10)*Sigma[2, 2]), n = 101,
#'                xvals = NULL, yvals = NULL, xlab = "x", ylab = "y",
#'                zlab = "f(x,y)", col = "gray", tol = 1e-6, ...)
#' @param mu length 2 numeric vector.
#' @param Sigma symmetric, positive-definite square matrix of order 2.
#' @param beta numeric. If \code{distribution = "mggd"}, the shape parameter of the MGGD.
#' \code{NULL} if \code{dist} is \code{"mcd"} or \code{"mtd"}.
#' @param nu numeric. If \code{distribution = "mtd"}, the degrees of freedom of the MTD.
#' \code{NULL} if \code{distribution} is \code{"mggd"} or \code{"mcd"}.
#' @param distribution the probability distribution. It can be \code{"mggd"} (multivariate generalized Gaussian distribution) \code{"mcd"} (multivariate Cauchy) or \code{"mtd"} (multivariate \eqn{t}).
#' @param xlim,ylim x-and y- limits.
#' @param n A one or two element vector giving the number of steps in the x and y grid, passed to \code{\link[rgl]{plot3d.function}}.
#' @param xvals,yvals The values at which to evaluate \code{x} and \code{y}. If used, \code{xlim} and/or \code{ylim} are ignored.
#' @param xlab,ylab,zlab The axis labels.
#' @param col The color to use for the plot. See \code{\link[rgl]{plot3d.function}}.
#' @param tol tolerance (relative to largest variance) for numerical lack of positive-definiteness in Sigma, for the estimation of the density. See \code{\link{dmggd}}, \code{\link{dmcd}} or \code{\link{dmtd}}.
#' @param ... Additional arguments to pass to \code{\link[rgl]{plot3d.function}}.
#' @return Returns invisibly the probability density function.
#'
#' @author Pierre Santagostini, Nizar Bouhlel
#' @references E. Gomez, M. Gomez-Villegas, H. Marin. A Multivariate Generalization of the Power Exponential Family of Distribution.
#' Commun. Statist. 1998, Theory Methods, col. 27, no. 23, p 589-600.
#' \doi{10.1080/03610929808832115}
#' 
#' S. Kotz and Saralees Nadarajah (2004), Multivariate \eqn{t} Distributions and Their Applications, Cambridge University Press.
#'
#' @seealso \code{\link{contourmvd}}: contour plot of a bivariate generalised Gaussian, Cauchy or \eqn{t} density.
#' 
#' \code{\link{dmggd}}: Probability density of a multivariate generalised Gaussian distribution.
#' 
#' \code{\link{dmcd}}: Probability density of a multivariate Cauchy distribution.
#' 
#' \code{\link{dmtd}}: Probability density of a multivariate \eqn{t} distribution.
#'
#' @examples
#' mu <- c(1, 4)
#' Sigma <- matrix(c(0.8, 0.2, 0.2, 0.2), nrow = 2)
#' 
#' # Bivariate generalised Gaussian distribution
#' beta <- 0.74
#' plotmvd(mu, Sigma, beta = beta, distribution = "mggd")
#' 
#' \donttest{
#' # Bivariate Cauchy distribution
#' plotmvd(mu, Sigma, distribution = "mcd")
#' 
#' # Bivariate t distribution
#' nu <- 2
#' plotmvd(mu, Sigma, nu = nu, distribution = "mtd")
#' }
#' 
#' @import rgl
#' @importFrom rgl plot3d
#' @export

plotmvd <- function(mu, Sigma, beta = NULL, nu = NULL,
                    distribution = c("mggd", "mcd", "mtd"),
                    xlim = c(mu[1] + c(-10, 10)*Sigma[1, 1]),
                    ylim = c(mu[2] + c(-10, 10)*Sigma[2, 2]), n = 101,
                    xvals = NULL, yvals = NULL, xlab = "x", ylab = "y",
                    zlab = "f(x,y)", col = "gray", tol = 1e-6, ...) {
  
  if (length(mu)!=2 | nrow(Sigma)!=2 | ncol(Sigma)!=2)
    stop(paste("plotmvd only allows plotting a bivariate density.",
               "mu must be a length 2 numeric vector and Sigma must be a 2*2 square matrix.", sep = "\n"))
  
  distribution <- match.arg(distribution)
  
  if (distribution == "mggd" & is.null(beta))
    stop('If distribution="mggd", beta must be provided.')
  if (distribution == "mtd" & is.null(nu))
    stop('If distribution="mtd", nu must be provided.')
  
  # Estimation of the density
  f <- switch(distribution,
              "mggd" = function(x) dmggd(x, mu = mu, Sigma = Sigma, beta = beta, tol = tol),
              "mcd" = function(x) dmcd(x, mu = mu, Sigma = Sigma, tol = tol),
              "mtd" = function(x) dmtd(x, nu = nu, mu = mu, Sigma = Sigma, tol = tol)
  )
  ff <- function(x, y) sapply(1:length(x), function(i) as.numeric(f(c(x[i], y[i]))))
  
  if (length(n) == 1)
    n <- rep(n, 2)
  if (is.null(xvals))
    xvals = seq.int(min(xlim), max(xlim), length.out = n[1])
  if (is.null(yvals))
    yvals = seq.int(min(ylim), max(ylim), length.out = n[2])
  
  # Plot
  plot3d(ff, xlim = xlim, ylim = ylim, n = n, xvals = xvals, yvals = yvals,
         xlab = xlab, ylab = ylab, zlab = zlab, col = col, ...)
  
  return(invisible(f))
}
