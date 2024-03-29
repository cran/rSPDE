% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inlabru_rspde.R
\name{bru_get_mapper.inla_rspde}
\alias{bru_get_mapper.inla_rspde}
\alias{ibm_n.bru_mapper_inla_rspde}
\alias{ibm_values.bru_mapper_inla_rspde}
\alias{ibm_jacobian.bru_mapper_inla_rspde}
\title{rSPDE inlabru mapper}
\usage{
bru_get_mapper.inla_rspde(model, ...)

ibm_n.bru_mapper_inla_rspde(mapper, ...)

ibm_values.bru_mapper_inla_rspde(mapper, ...)

ibm_jacobian.bru_mapper_inla_rspde(mapper, input, ...)
}
\arguments{
\item{model}{An \code{inla_rspde} for which to construct or extract a mapper}

\item{\dots}{Arguments passed on to other methods}

\item{mapper}{A \code{bru_mapper_inla_rspde} object}

\item{input}{The values for which to produce a mapping matrix}
}
\description{
rSPDE inlabru mapper
}
\examples{
\donttest{ #devel version
if (requireNamespace("INLA", quietly = TRUE) && 
     requireNamespace("inlabru", quietly = TRUE)){
library(INLA)
library(inlabru)

set.seed(123)
m <- 100
loc_2d_mesh <- matrix(runif(m * 2), m, 2)
mesh_2d <- inla.mesh.2d(
  loc = loc_2d_mesh,
  cutoff = 0.05,
  max.edge = c(0.1, 0.5)
)
sigma <- 1
range <- 0.2
nu <- 0.8
kappa <- sqrt(8 * nu) / range
op <- matern.operators(
  mesh = mesh_2d, nu = nu,
  range = range, sigma = sigma, m = 2,
  parameterization = "matern"
)
u <- simulate(op)
A <- inla.spde.make.A(
  mesh = mesh_2d,
  loc = loc_2d_mesh
)
sigma.e <- 0.1
y <- A \%*\% u + rnorm(m) * sigma.e
y <- as.vector(y)

data_df <- data.frame(y=y, x1 = loc_2d_mesh[,1],
                       x2 = loc_2d_mesh[,2])
coordinates(data_df) <- c("x1", "x2")
rspde_model <- rspde.matern(
  mesh = mesh_2d,
  nu_upper_bound = 2
)

cmp <- y ~ Intercept(1) + 
           field(coordinates, model = rspde_model)


rspde_fit <- bru(cmp, data = data_df)
summary(rspde_fit)
}
#devel.tag
}
}
