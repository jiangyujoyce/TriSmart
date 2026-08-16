\name{classify_missingness_zones}
\alias{classify_missingness_zones}
\title{Classify Features into TriSmart Missingness Zones}
\description{
Combines a smooth abundance/missingness trend with a three-component Gaussian
mixture model to assign each feature to one of three TriSmart zones.
}
\usage{
classify_missingness_zones(
  X,
  protein_names = NULL,
  gmm_model = "VEE",
  residual_quantile = 0.9,
  seed = 123
)
}
\arguments{
\item{X}{Numeric matrix or all-numeric data frame with samples in rows and at
least 15 named features in columns.}
\item{protein_names}{Optional unique feature labels. Defaults to column names.}
\item{gmm_model}{Covariance model name passed to \code{mclust::Mclust()}.}
\item{residual_quantile}{Quantile in \code{(0, 1)} for residual refinement.}
\item{seed}{Finite numeric seed.}
}
\value{
A list with \code{protein_zone_table}, fitted GAM and GMM objects,
\code{cluster_summary}, \code{residual_threshold}, and \code{lambda}.
}
\examples{
data(exampledata)
keep <- unique(round(seq(1, ncol(exampledata), length.out = 120)))
x <- exampledata[, keep, drop = FALSE]
\donttest{
if (requireNamespace("mclust", quietly = TRUE)) {
  zones <- classify_missingness_zones(x, seed = 2026)
  table(zones$protein_zone_table$final_zone)
}
}
}
