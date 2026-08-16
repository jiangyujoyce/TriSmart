\name{impute_mnar_by_col}
\alias{impute_mnar_by_col}
\title{Impute Low-Abundance Features}
\description{
Fills missing values only in selected low-abundance columns using a
left-censoring method while preserving observed and non-selected values.
}
\usage{
impute_mnar_by_col(
  X_imp,
  low_cols,
  method = c("QRILC", "MinProb", "Down-Shift", "MinDet", "GSimp"),
  q = 0.01,
  tune_sigma = 1,
  gsimp_iters_each = 50L,
  gsimp_iters_all = 10L,
  gsimp_cores = 1L,
  imputelcmd_axis = c("proteins", "samples"),
  seed = 12345
)
}
\arguments{
\item{X_imp}{Numeric sample-by-feature matrix.}
\item{low_cols}{Character vector naming the low-abundance columns to fill.}
\item{method}{One of \code{"QRILC"}, \code{"MinProb"},
\code{"MinDet"}, \code{"Down-Shift"}, or \code{"GSimp"}.}
\item{q}{Lower quantile in \code{(0, 1)} used by MinDet and MinProb.}
\item{tune_sigma}{Non-negative scale multiplier for QRILC or MinProb.}
\item{gsimp_iters_each, gsimp_iters_all, gsimp_cores}{Controls passed to
\code{GS_impute()}.}
\item{imputelcmd_axis}{Feature orientation. The recommended default is
\code{"proteins"}.}
\item{seed}{Finite numeric seed.}
}
\value{A numeric matrix with selected missing entries filled.}
\examples{
x <- matrix(
  c(10, 8, NA, 6, 4, 3, 2, 1),
  nrow = 4,
  dimnames = list(paste0("sample", 1:4), c("low", "reference"))
)
out <- impute_mnar_by_col(x, "low", method = "MinDet", seed = 1)
out
}
