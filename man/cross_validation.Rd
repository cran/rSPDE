% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inlabru_rspde.R
\name{cross_validation}
\alias{cross_validation}
\title{Perform cross-validation on a list of fitted models.}
\usage{
cross_validation(
  models,
  model_names = NULL,
  scores = c("mae", "mse", "crps", "scrps", "dss"),
  cv_type = c("k-fold", "loo", "lpo"),
  weight_thr = NULL,
  k = 5,
  percentage = 20,
  number_folds = 10,
  n_samples = 1000,
  return_scores_folds = FALSE,
  orientation_results = c("negative", "positive"),
  include_best = TRUE,
  train_test_indexes = NULL,
  return_train_test = FALSE,
  return_post_samples = FALSE,
  return_true_test_values = FALSE,
  parallelize_RP = FALSE,
  n_cores_RP = parallel::detectCores() - 1,
  true_CV = TRUE,
  save_settings = FALSE,
  print = TRUE,
  fit_verbose = FALSE
)
}
\arguments{
\item{models}{A fitted model obtained from calling the \code{bru()} function or a list of fitted models.
All models in the list must have the same number of likelihoods and must be fitted to
identical datasets.}

\item{model_names}{A vector containing the names of the models to appear in the returned \code{data.frame}. If \code{NULL}, the names will be of the form \verb{Model 1}, \verb{Model 2}, and so on. By default, it will try to obtain the name from the models list.}

\item{scores}{A vector containing the scores to be computed. The options are "mse", "crps", "scrps", "dss", "wcrps" and "swcrps". By default, all scores are computed.}

\item{cv_type}{The type of the folding to be carried out. The options are \code{k-fold} for \code{k}-fold cross-validation, in which case the parameter \code{k} should be provided,
\code{loo}, for leave-one-out and \code{lpo} for leave-percentage-out, in this case, the parameter \code{percentage} should be given, and also the \code{number_folds}
with the number of folds to be done. The default is \code{k-fold}.}

\item{weight_thr}{When computing "wcrps" or "swcrps", the threshold to be used to compute the weights. Must be supplied if any of these scores are requested. No default value is provided.}

\item{k}{The number of folds to be used in \code{k}-fold cross-validation. Will only be used if \code{cv_type} is \code{k-fold}.}

\item{percentage}{The percentage (from 1 to 99) of the data to be used to train the model. Will only be used if \code{cv_type} is \code{lpo}.}

\item{number_folds}{Number of folds to be done if \code{cv_type} is \code{lpo}.}

\item{n_samples}{Number of samples to compute the posterior statistics to be used to compute the scores.}

\item{return_scores_folds}{If \code{TRUE}, the scores for each fold will also be returned.}

\item{orientation_results}{character vector. The options are "negative" and "positive". If "negative", the smaller the scores the better. If "positive", the larger the scores the better.}

\item{include_best}{Should a row indicating which model was the best for each score be included?}

\item{train_test_indexes}{A list where each element corresponds to a fold. Each fold contains:
\itemize{
\item \code{train}: A list of training index vectors, one for each likelihood.
\item \code{test}: A list of test index vectors, one for each likelihood, with the same length as \code{train}.
This list is typically obtained by setting the argument \code{return_train_test} to \code{TRUE}.
When supplying \code{train_test_indexes}, the \code{cv_type}, \code{k}, \code{percentage} and \code{number_folds} arguments are ignored.
}}

\item{return_train_test}{Logical. Should the training and test indexes be returned? If 'TRUE' the train and test indexes will the 'train_test' element of the returned list.}

\item{return_post_samples}{If \code{TRUE} the posterior samples will be included in the returned list.}

\item{return_true_test_values}{If \code{TRUE} the true test values will be included in the returned list.}

\item{parallelize_RP}{Logical. Should the computation of CRPS and SCRPS (and for some cases, DSS) be parallelized?}

\item{n_cores_RP}{Number of cores to be used if \code{parallelize_rp} is \code{TRUE}.}

\item{true_CV}{Should a \code{TRUE} cross-validation be performed? If \code{TRUE} the models will be fitted on the training dataset. If \code{FALSE}, the parameters will be kept fixed at the ones obtained in the result object.}

\item{save_settings}{Logical. If \code{TRUE}, the settings used in the cross-validation will also be returned.}

\item{print}{Should partial results be printed throughout the computation?}

\item{fit_verbose}{Should INLA's run during cross-validation be verbose?}
}
\value{
A data.frame with the fitted models and the corresponding scores.
}
\description{
Obtain several scores for a list of fitted models according
to a folding scheme.
}
