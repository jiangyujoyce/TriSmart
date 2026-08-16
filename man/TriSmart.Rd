\name{TriSmart}
\alias{TriSmart}
\alias{tri_smart}
\title{Run the TriSmart Hybrid Imputation Pipeline}
\description{
Classifies features into three abundance/missingness zones and applies a
three-stage imputation workflow while preserving every originally observed
value.
}
\usage{
TriSmart(
  X_miss,
  firstImpute,
  secondImpute,
  workers = NULL,
  max_workers = 8L,
  seed = 1L,
  residual_quantile = 0.9,
  gmm_model = "VEE",
  min_observed = 10L,
  max_predictors = NULL,
  glmnet_alpha = 0.5,
  glmnet_nlambda = 100L,
  knn_axis = c("proteins", "samples"),
  imputelcmd_axis = c("proteins", "samples"),
  rf_ntree = 100L,
  rf_maxiter = 10L,
  bpca_npcs = 5L,
  mnar_q = 0.01,
  mnar_tune_sigma = 1,
  gsimp_cores = 1L,
  verbose = FALSE,
  return_details = FALSE
)

tri_smart(
  X_miss,
  firstImpute,
  secondImpute,
  workers = NULL,
  max_workers = 8L,
  seed = 1L,
  residual_quantile = 0.9,
  gmm_model = "VEE",
  min_observed = 10L,
  max_predictors = NULL,
  glmnet_alpha = 0.5,
  glmnet_nlambda = 100L,
  knn_axis = c("proteins", "samples"),
  imputelcmd_axis = c("proteins", "samples"),
  rf_ntree = 100L,
  rf_maxiter = 10L,
  bpca_npcs = 5L,
  mnar_q = 0.01,
  mnar_tune_sigma = 1,
  gsimp_cores = 1L,
  verbose = FALSE,
  return_details = FALSE
)
}
\arguments{
\item{X_miss}{Numeric matrix or all-numeric data frame with samples in rows
and uniquely named proteins/features in columns.}
\item{firstImpute}{First-stage method: \code{"KNN"},
\code{"RandomForest"}, \code{"BPCA"}, or \code{"Ensemble"}.}
\item{secondImpute}{Third-stage method: \code{"QRILC"},
\code{"MinProb"}, \code{"MinDet"}, \code{"Down-Shift"}, or
\code{"GSimp"}. \code{"Down-Shift"} is a compatibility alias for
\code{"MinDet"}.}
\item{workers}{Number of workers for mixture-model targets. \code{NULL}
selects a conservative value automatically.}
\item{max_workers}{Upper bound on automatically or explicitly selected workers.}
\item{seed}{Finite numeric seed used reproducibly across all stages.}
\item{residual_quantile}{Quantile in \code{(0, 1)} used to identify large
deviations from the abundance/missingness trend.}
\item{gmm_model}{Covariance model name passed to \code{mclust::Mclust()}.}
\item{min_observed}{Minimum observed target values required for elastic-net
modeling.}
\item{max_predictors}{Maximum predictors per mixture-zone target, or
\code{NULL} for an automatic choice.}
\item{glmnet_alpha}{Elastic-net mixing parameter between zero and one.}
\item{glmnet_nlambda}{Number of lambda values for cross-validation.}
\item{knn_axis}{Orientation used by KNN. The default follows the package's
protein-column convention; \code{"samples"} is retained for compatibility.}
\item{imputelcmd_axis}{Whether imputeLCMD-style methods treat proteins or
samples as features. The recommended default is \code{"proteins"}.}
\item{rf_ntree}{Number of random-forest trees.}
\item{rf_maxiter}{Maximum missForest iterations.}
\item{bpca_npcs}{Requested number of BPCA components.}
\item{mnar_q}{Lower quantile used by MinDet and MinProb.}
\item{mnar_tune_sigma}{Non-negative scale multiplier for QRILC or MinProb.}
\item{gsimp_cores}{Number of cores passed to \code{GS_impute()}.}
\item{verbose}{Whether to report progress.}
\item{return_details}{If \code{TRUE}, return imputed data, zone assignments,
diagnostics, and settings; otherwise return only the imputed matrix.}
}
\value{
By default, a completely imputed numeric matrix with the same dimensions and
dimnames as \code{X_miss}. If \code{return_details = TRUE}, a list containing
\code{imputed_data}, \code{zone_result}, \code{zone_columns},
\code{mixture_diagnostics}, and \code{settings}.
}
\details{
Zone 1 contains relatively high-abundance features and Zone 3 relatively
low-abundance features. Zone 2 contains mixed patterns. The chosen first-stage
method fills Zone 1 and Zone 2, elastic-net models refine eligible Zone 2
values, and the selected left-censoring method fills Zone 3.
}
\examples{
data(exampledata)
keep <- unique(round(seq(1, ncol(exampledata), length.out = 120)))
x <- exampledata[, keep, drop = FALSE]

\donttest{
if (requireNamespace("mclust", quietly = TRUE) &&
    requireNamespace("impute", quietly = TRUE) &&
    requireNamespace("glmnet", quietly = TRUE)) {
  fit <- TriSmart(
    x,
    firstImpute = "KNN",
    secondImpute = "MinDet",
    workers = 1,
    seed = 2026,
    return_details = TRUE
  )
  anyNA(fit$imputed_data)
  table(fit$zone_result$protein_zone_table$final_zone)
}
}
}
