% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spacetime.R
\name{update.spacetimeobj}
\alias{update.spacetimeobj}
\title{Update parameters of spacetimeobj objects}
\usage{
\method{update}{spacetimeobj}(object, kappa = NULL, sigma = NULL, gamma = NULL, rho = NULL, ...)
}
\arguments{
\item{object}{Space-time object created by \code{\link[=spacetime.operators]{spacetime.operators()}}}

\item{kappa}{kappa value to be updated.}

\item{sigma}{sigma value to be updated.}

\item{gamma}{gamma value to be updated.}

\item{rho}{rho value to be updated.}

\item{...}{currently not used.}
}
\value{
An object of type spacetimeobj with updated parameters.
}
\description{
Function to change the parameters of a spacetimeobj object
}
\examples{
s <- seq(from = 0, to = 20, length.out = 101)
t <- seq(from = 0, to = 20, length.out = 31)

op_cov <- spacetime.operators(space_loc = s, time_loc = t,
                             kappa = 5, sigma = 10, alpha = 1,
                             beta = 2, rho = 1, gamma = 0.05)
op_cov <- update(op_cov, kappa = 4, 
                             sigma = 2, gamma = 0.1)                              
}
