% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inla_rspde_1d.R
\name{bru_get_mapper.inla_rspde_matern1d}
\alias{bru_get_mapper.inla_rspde_matern1d}
\alias{ibm_n.bru_mapper_inla_rspde_matern1d}
\alias{ibm_values.bru_mapper_inla_rspde_matern1d}
\alias{ibm_jacobian.bru_mapper_inla_rspde_matern1d}
\title{rSPDE stationary inlabru mapper}
\usage{
bru_get_mapper.inla_rspde_matern1d(model, ...)

ibm_n.bru_mapper_inla_rspde_matern1d(mapper, ...)

ibm_values.bru_mapper_inla_rspde_matern1d(mapper, ...)

ibm_jacobian.bru_mapper_inla_rspde_matern1d(mapper, input, ...)
}
\arguments{
\item{model}{An \code{inla_rspde_matern1d} object for which to construct or extract a mapper}

\item{\dots}{Arguments passed on to other methods}

\item{mapper}{A \code{bru_mapper_inla_rspde_matern1d} object}

\item{input}{The values for which to produce a mapping matrix}
}
\description{
rSPDE stationary inlabru mapper
}
