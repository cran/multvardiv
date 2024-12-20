contourmvd <- function(mu, Sigma, beta = NULL, nu = NULL,
                        distribution = c("mggd", "mcd", "mtd"),
                        xlim = c(mu[1] + c(-10, 10)*Sigma[1, 1]),
                        ylim = c(mu[2] + c(-10, 10)*Sigma[2, 2]),
                        zlim = NULL, npt = 30, nx = npt, ny = npt,
                        main = NULL, sub = NULL, nlevels = 10,
                        levels = pretty(zlim, nlevels),
                        tol = 1e-6, ...) {
  #' Contour Plot of a Bivariate Density
  #' 
  #' Contour plot of the probability density of a multivariate distribution with 2 variables:
  #' \itemize{
  #' \item generalized Gaussian distribution (MGGD) with mean vector \code{mu}, dispersion matrix \code{Sigma} and shape parameter \code{beta}
  #' \item Cauchy distribution (MCD) with location parameter \code{mu} and scatter matrix \code{Sigma}
  #' \item \eqn{t} distribution (MTD) with location parameter \code{mu}, scatter matrix \code{Sigma} and degrees of freedom \code{nu}
  #' }
  #' This function uses the \code{\link{contour}} function.
  #'
  #' @aliases contourmvd
  #'
  #' @usage contourmvd(mu, Sigma, beta = NULL, nu = NULL,
  #'                   distribution = c("mggd", "mcd", "mtd"),
  #'                   xlim = c(mu[1] + c(-10, 10)*Sigma[1, 1]),
  #'                   ylim = c(mu[2] + c(-10, 10)*Sigma[2, 2]),
  #'                   zlim = NULL, npt = 30, nx = npt, ny = npt,
  #'                   main = NULL, sub = NULL, nlevels = 10,
  #'                   levels = pretty(zlim, nlevels), tol = 1e-6, ...)
  #' @param mu length 2 numeric vector.
  #' @param Sigma symmetric, positive-definite square matrix of order 2. The dispersion matrix.
  #' @param beta numeric. If \code{distribution = "mggd"}, the shape parameter of the MGGD.
  #' \code{NULL} if \code{dist} is \code{"mcd"} or \code{"mtd"}.
  #' @param nu numeric. If \code{distribution = "mtd"}, the degrees of freedom of the MTD.
  #' \code{NULL} if \code{distribution} is \code{"mggd"} or \code{"mcd"}.
  #' @param distribution character string. The probability distribution. It can be \code{"mggd"} (multivariate generalized Gaussian distribution) \code{"mcd"} (multivariate Cauchy) or \code{"mtd"} (multivariate \eqn{t}).
  #' @param main,sub main and sub title, as for \code{\link{title}}.
  #' If omitted, the main title is set to `"Multivariate generalised Gaussian density"`,
  #' `"Multivariate Cauchy density"` or `"Multivariate t density"`.
  #' @param xlim,ylim x-and y- limits.
  #' @param zlim z- limits. If NULL, it is the range of the values of the density on the x and y values within `xlim` and `ylim`.
  #' @param npt number of points for the discretisation.
  #' @param nx,ny number of points for the discretisation among the x- and y- axes.
  #' @param nlevels,levels arguments to be passed to the \code{\link{contour}} function.
  #' @param tol tolerance (relative to largest variance) for numerical lack of positive-definiteness in Sigma, for the estimation of the density. See \code{\link{dmggd}}, \code{\link{dmcd}} or \code{\link{dmtd}}.
  #' @param ... additional arguments to \code{\link{plot.window}}, \code{\link{title}}, \code{\link{Axis}} and \code{\link{box}}, typically \link{graphical parameters} such as \code{cex.axis}.
  #' @return Returns invisibly the probability density function.
  #'
  #' @author Pierre Santagostini, Nizar Bouhlel
  #' @references E. Gomez, M. Gomez-Villegas, H. Marin. A Multivariate Generalization of the Power Exponential Family of Distribution.
  #' Commun. Statist. 1998, Theory Methods, col. 27, no. 23, p 589-600.
  #' \doi{10.1080/03610929808832115}
  #' 
  #' S. Kotz and Saralees Nadarajah (2004), Multivariate \eqn{t} Distributions and Their Applications, Cambridge University Press.
  #'
  #' @seealso \code{\link{plotmvd}}: plot of a bivariate generalised Gaussian, Cauchy or \eqn{t} density.
  #' 
  #' \code{\link{dmggd}}: probability density of a multivariate generalised Gaussian distribution.
  #' 
  #' \code{\link{dmcd}}: probability density of a multivariate Cauchy distribution.
  #' 
  #' \code{\link{dmtd}}: probability density of a multivariate \eqn{t} distribution.
  #'
  #' @examples
  #' mu <- c(1, 4)
  #' Sigma <- matrix(c(0.8, 0.2, 0.2, 0.2), nrow = 2)
  #' 
  #' # Bivariate generalized Gaussian distribution
  #' beta <- 0.74
  #' contourmvd(mu, Sigma, beta = beta, distribution = "mggd")
  #' 
  #' # Bivariate Cauchy distribution
  #' contourmvd(mu, Sigma, distribution = "mcd")
  #' 
  #' # Bivariate t distribution
  #' nu <- 1
  #' contourmvd(mu, Sigma, nu = nu, distribution = "mtd")
  #'
  #' @importFrom graphics contour
  #' @importFrom graphics par
  #' @export
  
  if (length(mu)!=2 | nrow(Sigma)!=2 | ncol(Sigma)!=2)
    stop(paste("contourmggd only allows plotting a generalised Gaussian density with 2 variables.",
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
  
  x <- seq(xlim[1], xlim[2], length = nx)
  y <- seq(ylim[1], ylim[2], length = ny)
  z <- outer(x, y, ff)
  if (is.null(zlim)) zlim <- range(z)
  
  if (is.null(main)) {
    main <- c(mggd = "generalised Gaussian", mcd = "Cauchy", mtd = "t")[distribution]
    main <- paste("Multivariate", main, "density")
  }
  
  # Plot
  contour(x, y, z, nlevels = nlevels, levels = levels, labels = NULL,
          xlim = xlim, ylim = ylim, zlim = zlim,  labcex = 0.6,
          drawlabels = TRUE, method = "flattest", vfont = NULL, axes = TRUE,
          frame.plot = TRUE, col = par("fg"), lty = par("lty"),
          lwd = par("lwd"), add = FALSE, main = main, ...)
  
  return(invisible(f))
}
