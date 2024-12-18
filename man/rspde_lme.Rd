% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rspde_lme.R
\name{rspde_lme}
\alias{rspde_lme}
\title{rSPDE linear mixed effects models}
\usage{
rspde_lme(
  formula,
  loc,
  loc_time = NULL,
  data,
  model = NULL,
  repl = NULL,
  which_repl = NULL,
  optim_method = "L-BFGS-B",
  use_data_from_graph = TRUE,
  starting_values_latent = NULL,
  start_sigma_e = NULL,
  start_alpha = NULL,
  alpha = NULL,
  start_nu = NULL,
  nu = NULL,
  nu_upper_bound = 4,
  rspde_order = NULL,
  parallel = FALSE,
  n_cores = parallel::detectCores() - 1,
  optim_controls = list(),
  improve_hessian = FALSE,
  hessian_args = list()
)
}
\arguments{
\item{formula}{Formula object describing the relation between the response
variables and the fixed effects. If the response variable is a matrix, each
column of the matrix will be treated as a replicate.}

\item{loc}{A vector with the names of the columns in \code{data} that contain the
observation locations, or a \code{matrix} or a \code{data.frame} containing the
observation locations. If the model is of class \code{metric_graph}, the locations
must be either a \code{matrix} or a \code{data.frame} with two columns, or a character
vector with the names of the two columns. The first column being the number of
the edge, and the second column being the normalized position on the edge.
If the model is a 2d model, \code{loc} must be either a \code{matrix} or \code{data.frame}
with two columns or a character vector with the name of the two columns that
contain the location, the first entry corresponding to the \code{x} entry and the
second corresponding to the \code{y} entry.}

\item{loc_time}{For spatio-temporal models, the name of the column in \code{data} that
is the time variable, or a \code{matrix} or \code{vector} containing the observation time
points.}

\item{data}{A \code{data.frame} containing the data to be used.}

\item{model}{Object generated by \code{matern.operators()}, \code{spde.matern.operators()}
or \code{spacetime.operators()}. If \code{NULL}, simple linear regression will be performed.}

\item{repl}{Vector indicating the replicate of each observation.
If \code{NULL} it will assume there is only one replicate. If the model is generated from graphs from
\code{metric_graph} class and \code{use_data_from_graph} is \code{TRUE}, \code{repl} needs to be the name of the
column inside the metric graph data that contains the replicate. If \code{NULL} it will assume there is only
one replicate.}

\item{which_repl}{Which replicates to use? If \code{NULL} all replicates will be used.}

\item{optim_method}{The method to be used with \code{optim} function.}

\item{use_data_from_graph}{Logical. Only for models generated from graphs from
\code{metric_graph} class. In this case, should the data, the locations and the
replicates be obtained from the graph object?}

\item{starting_values_latent}{A vector containing the starting values for the
latent model. If the latent model is generated by \code{matern.operators()}, then
the vector should be on the form c(tau,kappa). If the model is generated by
\code{spde.matern.operators()}, the vector should contain the nonstationary parameters.
If the model is generated by \code{spacetime.operators()}, the vector should be on
the form c(kappa,sigma,gamma,rho).}

\item{start_sigma_e}{Starting value for the standard deviation of the
measurement error.}

\item{start_alpha}{Starting value for the smoothness parameter of spatial
models. Will be used if \code{start_nu} is not given. Not used for spatio-temporal
models.}

\item{alpha}{If \code{NULL}, the smoothness parameter will be estimated for spatial
models, otherwise it is kept fixed at the provided value. Will be used if \code{nu}
is not given. Not used for spatio-temporal models.
returned as component of the returned value.}

\item{start_nu}{Starting value for the smoothness parameter of spatial models.
Not used for spatio-temporal models.}

\item{nu}{If \code{NULL}, the smoothness parameter will be estimated for spatial
models, otherwise the smoothness parameter will be kept fixed at the provided
value. Not used for spatio-temporal models.}

\item{nu_upper_bound}{A parameter that limits the maximum value that nu can
assume. Not used for spatio-temporal models.}

\item{rspde_order}{The order of the rational approximation to be used while
fitting the model. If not given, the order from the model object will be used.
Not used for spatio-temporal models.}

\item{parallel}{logical. Indicating whether to use optimParallel or not.}

\item{n_cores}{Number of cores to be used if parallel is true.}

\item{optim_controls}{Additional controls to be passed to \code{optim} or \code{optimParallel}.}

\item{improve_hessian}{Should a more precise estimate of the hessian be obtained?
Turning on might increase the overall time.}

\item{hessian_args}{List of controls to be used if \code{improve_hessian} is \code{TRUE}.
The list can contain the arguments to be passed to the \code{method.args} argument
in the \code{numDeriv::hessian} function. See the help of the \code{hessian} function in
\code{numDeriv} package for details. Observe that it only accepts the "Richardson"
method for now, the method "complex" is not supported.}
}
\value{
A list containing the fitted model.
}
\description{
Fitting linear mixed effects model with latent Whittle-Matern models.
}
