\name{exampledata}
\alias{exampledata}
\docType{data}
\title{Example Quantitative Proteomics Data with Missing Values}
\description{
A numeric matrix supplied with TriSmart for illustrating the hybrid
imputation workflow. Rows represent samples and columns represent proteins or
features. Missing measurements are stored as \code{NA}.
}
\usage{data(exampledata)}
\format{
A numeric matrix with 20 rows (samples) and 6,611 columns (proteins/features).
Row names identify samples and column names identify features.
}
\source{Example dataset supplied with the original TriSmart analysis code.}
\examples{
data(exampledata)
dim(exampledata)
mean(is.na(exampledata))

# A smaller matrix suitable for an introductory analysis.
keep <- unique(round(seq(1, ncol(exampledata), length.out = 120)))
example_subset <- exampledata[, keep, drop = FALSE]
dim(example_subset)
}
\keyword{datasets}
