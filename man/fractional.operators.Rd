% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fractional.operators.R
\name{fractional.operators}
\alias{fractional.operators}
\title{Rational approximations of fractional operators}
\usage{
fractional.operators(L, beta, C, scale.factor, m = 1, tau = 1)
}
\arguments{
\item{L}{A finite element discretization of the operator \eqn{L}{L}.}

\item{beta}{The positive fractional power.}

\item{C}{The mass matrix of the finite element discretization.}

\item{scale.factor}{A constant \eqn{c} is a lower bound for the the smallest
eigenvalue of the non-discretized operator \eqn{L}{L}.}

\item{m}{The order of the rational approximation, which needs to be a
positive integer. The default value is 1.
Higer values gives a more accurate approximation, which are more
computationally expensive to use for inference. Currently, the largest value
of m that is implemented is 4.}

\item{tau}{The constant or vector that scales the variance of the solution.
The default value is 1.}
}
\value{
\code{fractional.operators} returns an object of class "rSPDEobj".
This object contains the following quantities:
\item{Pl}{The operator \eqn{P_l}.}
\item{Pr}{The operator \eqn{P_r}.}
\item{C}{The mass lumped mass matrix.}
\item{Ci}{The inverse of \code{C}.}
\item{m}{The order of the rational approximation.}
\item{beta}{The fractional power.}
\item{type}{String indicating the type of approximation.}
\item{Q}{The matrix \code{t(Pl) \%*\% solve(C,Pl)}.}
\item{type}{String indicating the type of approximation.}
\item{Pl.factors}{List with elements that can be used to assemble \eqn{P_l}.}
\item{Pr.factors}{List with elements that can be used to assemble \eqn{P_r}.}
}
\description{
\code{fractional.operators} is used for computing an approximation,
which can be used for inference and simulation, of the fractional SPDE
\deqn{L^\beta (\tau u(s)) = W.}
Here \eqn{L} is a differential operator, \eqn{\beta>0} is
the fractional power, \eqn{\tau} is a positive scalar or vector that
scales the variance of the solution \eqn{u}, and \eqn{W} is white noise.
}
\details{
The approximation is based on a rational approximation of
the fractional operator, resulting in an
approximate model on the form \deqn{P_l u(s) = P_r W,}
where \eqn{P_j = p_j(L)} are non-fractional operators defined in terms of
polynomials \eqn{p_j} for \eqn{j=l,r}. The order of \eqn{p_r} is given by
\code{m} and the order of \eqn{p_l} is \eqn{m + m_\beta}
where \eqn{m_\beta} is the integer part of \eqn{\beta} if \eqn{\beta>1} and
\eqn{m_\beta = 1} otherwise.

The discrete approximation can be written as \eqn{u = P_r x} where
\eqn{x \sim N(0,Q^{-1})}{x ~ N(0,Q^{-1})} and \eqn{Q = P_l^T C^{-1} P_l}.
Note that the matrices \eqn{P_r} and \eqn{Q} may be be ill-conditioned
for \eqn{m>1}. In this case, the methods in \code{\link[=operator.operations]{operator.operations()}}
should be used for operations involving the matrices, since these methods
are more numerically stable.
}
\examples{
# Compute rational approximation of a Gaussian process with a
# Matern covariance function on R
kappa <- 10
sigma <- 1
nu <- 0.8

# create mass and stiffness matrices for a FEM discretization
x <- seq(from = 0, to = 1, length.out = 101)
fem <- rSPDE.fem1d(x)

# compute rational approximation of covariance function at 0.5
tau <- sqrt(gamma(nu) / (sigma^2 * kappa^(2 * nu) *
  (4 * pi)^(1 / 2) * gamma(nu + 1 / 2)))
op <- fractional.operators(
  L = fem$G + kappa^2 * fem$C, beta = (nu + 1 / 2) / 2,
  C = fem$C, scale.factor = kappa^2, tau = tau
)

v <- t(rSPDE.A1d(x, 0.5))
c.approx <- Sigma.mult(op, v)

# plot the result and compare with the true Matern covariance
plot(x, matern.covariance(abs(x - 0.5), kappa, nu, sigma),
  type = "l", ylab = "C(h)",
  xlab = "h", main = "Matern covariance and rational approximations"
)
lines(x, c.approx, col = 2)

}
\seealso{
\code{\link[=matern.operators]{matern.operators()}}, \code{\link[=spde.matern.operators]{spde.matern.operators()}},
\code{\link[=matern.operators]{matern.operators()}}
}
