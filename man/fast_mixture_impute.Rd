\name{fast_mixture_impute}
\alias{fast_mixture_impute}
\title{Refine Imputations for Mixed-Missingness Features}
\description{
Uses elastic-net prediction to refine first-stage values for eligible Zone 2
features, with the first-stage result as a safe fallback.
}
\usage{
fast_mixture_impute(
  X_miss,
  X_impute_step1,
  mixture_idx,
  high_idx,
  min_observed = 10L,
  max_predictors = NULL,
  alpha = 0.5,
  nlambda = 100L,
  workers = NULL,
  max_workers = 8L,
  parallel_min_tasks = 25L,
  seed = 1234,
  return_details = FALSE
)
}
\arguments{
\item{X_miss}{Original numeric sample-by-feature matrix.}
\item{X_impute_step1}{Full matrix or high/mixture subset after first-stage
imputation.}
\item{mixture_idx, high_idx}{Character vectors naming Zone 2 and Zone 1
columns.}
\item{min_observed}{Minimum observed target values required for modeling.}
\item{max_predictors}{Maximum selected predictors per target.}
\item{alpha}{Elastic-net mixing parameter in \code{[0, 1]}.}
\item{nlambda}{Number of lambda values for cross-validation.}
\item{workers, max_workers}{Requested and maximum parallel worker counts.}
\item{parallel_min_tasks}{Minimum eligible targets before parallel execution.}
\item{seed}{Finite numeric seed.}
\item{return_details}{Whether to include model diagnostics.}
}
\value{
An imputed matrix, or a list with \code{data} and \code{diagnostics} when
\code{return_details = TRUE}.
}
