\name{impute_high_abundance}
\alias{impute_high_abundance}
\title{Impute High- and Mixed-Abundance Features}
\description{
Performs the first TriSmart stage on selected Zone 1 and Zone 2 columns using
KNN, random forest, Bayesian PCA, or their mean ensemble.
}
\usage{
impute_high_abundance(
  X_miss,
  zone1_cols,
  zone2_cols = NULL,
  impmethod = c("KNN", "RandomForest", "BPCA", "Ensemble"),
  rf_ntree = 100L,
  rf_maxiter = 10L,
  bpca_npcs = 5L,
  knn_k = 10L,
  knn_rowmax = 0.9,
  knn_colmax = 0.9,
  knn_maxp = 1500L,
  knn_axis = c("proteins", "samples"),
  rf_workers = 1L,
  seed = 362436069,
  verbose = FALSE
)
}
\arguments{
\item{X_miss}{Numeric sample-by-feature matrix containing missing values.}
\item{zone1_cols, zone2_cols}{Character vectors of feature column names.}
\item{impmethod}{One of \code{"KNN"}, \code{"RandomForest"},
\code{"BPCA"}, or \code{"Ensemble"}.}
\item{rf_ntree, rf_maxiter}{Random-forest tuning controls.}
\item{bpca_npcs}{Requested number of BPCA components.}
\item{knn_k, knn_rowmax, knn_colmax, knn_maxp}{Arguments used by
\code{impute::impute.knn()}.}
\item{knn_axis}{KNN orientation; \code{"proteins"} is recommended.}
\item{rf_workers}{Validated worker count retained for interface consistency.}
\item{seed}{Finite numeric seed.}
\item{verbose}{Whether the selected backend should report progress.}
}
\value{A completely imputed matrix containing the selected columns only.}
\examples{
data(exampledata)
x <- exampledata[, 1:30, drop = FALSE]
\donttest{
if (requireNamespace("impute", quietly = TRUE)) {
  out <- impute_high_abundance(x, colnames(x), impmethod = "KNN")
  anyNA(out)
}
}
}
