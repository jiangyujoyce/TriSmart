library(TriSmart)

data(exampledata)
stopifnot(
  is.matrix(exampledata),
  is.numeric(exampledata),
  identical(dim(exampledata), c(20L, 6611L)),
  anyNA(exampledata),
  !is.null(rownames(exampledata)),
  !is.null(colnames(exampledata))
)

# A dependency-free check of the third-stage MinDet implementation.
x <- matrix(
  c(10, 8, NA, 6, 4, 3, 2, 1),
  nrow = 4,
  dimnames = list(paste0("sample", 1:4), c("low", "reference"))
)
observed <- !is.na(x)
out <- impute_mnar_by_col(x, low_cols = "low", method = "MinDet", seed = 1)
stopifnot(!anyNA(out), identical(out[observed], x[observed]))
