# Rebuild data/exampledata.rda from the original supplied workspace file.
source_file <- "Exampledata.Rdata"
source_environment <- new.env(parent = emptyenv())
loaded_objects <- load(source_file, envir = source_environment)

if (!identical(loaded_objects, "Mix")) {
  stop("Exampledata.Rdata must contain exactly one object named 'Mix'.")
}

exampledata <- source_environment$Mix
stopifnot(
  is.matrix(exampledata),
  is.numeric(exampledata),
  !is.null(rownames(exampledata)),
  !is.null(colnames(exampledata))
)

save(exampledata, file = "data/exampledata.rda", compress = "xz")
