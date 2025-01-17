% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rmtd.R
\name{rmtd}
\alias{rmtd}
\title{Simulate from a Multivariate \eqn{t} Distribution}
\usage{
rmtd(n, nu, mu, Sigma, tol = 1e-6)
}
\arguments{
\item{n}{integer. Number of observations.}

\item{nu}{numeric. The degrees of freedom.}

\item{mu}{length \eqn{p} numeric vector. The mean vector}

\item{Sigma}{symmetric, positive-definite square matrix of order \eqn{p}. The correlation matrix.}

\item{tol}{tolerance for numerical lack of positive-definiteness in Sigma (for \code{mvrnorm}, see Details).}
}
\value{
A matrix with \eqn{p} columns and \eqn{n} rows.
}
\description{
Produces one or more samples from the multivariate (\eqn{p} variables) \eqn{t} distribution (MTD)
with degrees of freedom \code{nu}, mean vector \code{mu} and
correlation matrix \code{Sigma}.
}
\details{
A sample from a MTD with parameters \eqn{\nu}, \eqn{\boldsymbol{\mu}} and \eqn{\Sigma}
can be generated using:
\deqn{\displaystyle{\mathbf{X} = \boldsymbol{\mu} + \mathbf{Y} \sqrt{\frac{\nu}{u}}}}
where \eqn{Y} is a random vector distributed among a centered Gaussian density
with covariance matrix \eqn{\Sigma} (generated using \code{\link[MASS]{mvrnorm}})
and \eqn{u} is distributed among a Chi-squared distribution with \eqn{\nu} degrees of freedom.
}
\examples{
nu <- 3
mu <- c(0, 1, 4)
Sigma <- matrix(c(1, 0.6, 0.2, 0.6, 1, 0.3, 0.2, 0.3, 1), nrow = 3)
x <- rmtd(10000, nu, mu, Sigma)
head(x)
dim(x)
mu; colMeans(x)
nu/(nu-2)*Sigma; var(x)

}
\references{
S. Kotz and Saralees Nadarajah (2004), Multivariate \eqn{t} Distributions and Their Applications, Cambridge University Press.
}
\seealso{
\code{\link{dmtd}}: probability density of a MTD.

\code{\link{estparmtd}}: estimation of the parameters of a MTD.
}
\author{
Pierre Santagostini, Nizar Bouhlel
}
