# TriSmart optimized, corrected, and portable implementation
#
# Input convention:
#   rows    = samples
#   columns = proteins/features
#
# Core packages are referenced with pkg::fun to avoid search-path conflicts.
# Optional first-stage packages are checked only when their methods are used.

# -----------------------------
# General helpers
# -----------------------------

.trismart_require <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(
      sprintf(
        "Package '%s' is required for this operation but is not installed.",
        package
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.trismart_seed_offset <- function(seed, offset = 0L) {
  seed <- as.double(seed)
  offset <- as.double(offset)
  if (length(seed) != 1L || is.na(seed) || !is.finite(seed)) {
    stop("seed must be a single finite number.", call. = FALSE)
  }
  if (length(offset) != 1L || is.na(offset) || !is.finite(offset)) {
    stop("offset must be a single finite number.", call. = FALSE)
  }
  value <- (seed + offset) %% .Machine$integer.max
  if (value <= 0) value <- value + 1
  as.integer(value)
}


.trismart_integer_scalar <- function(value, name, minimum = 1L) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || value != floor(value) || value < minimum ||
      value > .Machine$integer.max) {
    stop(
      sprintf("%s must be one whole number of at least %s.", name, minimum),
      call. = FALSE
    )
  }
  as.integer(value)
}

.trismart_with_seed <- function(seed, code) {
  if (is.null(seed)) {
    return(force(code))
  }

  seed <- .trismart_seed_offset(seed)
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }

  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(seed)
  force(code)
}

.trismart_validate_matrix <- function(X, min_rows = 1L, min_cols = 1L) {
  if (is.data.frame(X)) {
    numeric_columns <- vapply(X, is.numeric, logical(1))
    if (!all(numeric_columns)) {
      bad <- names(numeric_columns)[!numeric_columns]
      stop(
        sprintf(
          "All columns must be numeric. Non-numeric columns: %s",
          paste(bad, collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }

  X <- as.matrix(X)
  if (!is.numeric(X)) {
    stop("X must be a numeric matrix or an all-numeric data frame.", call. = FALSE)
  }
  storage.mode(X) <- "double"

  if (length(dim(X)) != 2L || nrow(X) < min_rows || ncol(X) < min_cols) {
    stop(
      sprintf(
        "X must contain at least %d samples (rows) and %d proteins (columns).",
        min_rows,
        min_cols
      ),
      call. = FALSE
    )
  }

  if (is.null(colnames(X))) {
    stop("X must have protein/feature column names.", call. = FALSE)
  }
  if (anyNA(colnames(X)) || any(trimws(colnames(X)) == "")) {
    stop("All protein/feature columns must have non-empty names.", call. = FALSE)
  }
  if (anyDuplicated(colnames(X))) {
    duplicates <- unique(colnames(X)[duplicated(colnames(X))])
    stop(
      sprintf(
        "Protein/feature names must be unique. Duplicates: %s",
        paste(duplicates, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (any(is.infinite(X), na.rm = TRUE)) {
    stop("X contains Inf or -Inf; replace or remove them before imputation.", call. = FALSE)
  }

  all_missing <- colSums(!is.na(X)) == 0L
  if (any(all_missing)) {
    stop(
      sprintf(
        paste0(
          "Completely missing proteins cannot be estimated without an external rule. ",
          "Remove them or supply values first: %s"
        ),
        paste(colnames(X)[all_missing], collapse = ", ")
      ),
      call. = FALSE
    )
  }

  all_missing_samples <- rowSums(!is.na(X)) == 0L
  if (any(all_missing_samples)) {
    sample_labels <- rownames(X)
    if (is.null(sample_labels)) sample_labels <- as.character(seq_len(nrow(X)))
    stop(
      sprintf(
        paste0(
          "Completely missing samples cannot be imputed reliably by this pipeline. ",
          "Remove them or supply values first: %s"
        ),
        paste(sample_labels[all_missing_samples], collapse = ", ")
      ),
      call. = FALSE
    )
  }

  X
}

.trismart_validate_columns <- function(columns, X, argument_name) {
  columns <- unique(as.character(columns))
  columns <- columns[!is.na(columns) & nzchar(columns)]
  unknown <- setdiff(columns, colnames(X))
  if (length(unknown)) {
    stop(
      sprintf(
        "%s contains unknown columns: %s",
        argument_name,
        paste(unknown, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  columns
}

.trismart_align_candidate <- function(original, candidate, method_name) {
  candidate <- as.matrix(candidate)
  if (!identical(dim(original), dim(candidate))) {
    stop(
      sprintf("%s returned an object with incompatible dimensions.", method_name),
      call. = FALSE
    )
  }

  original_rows <- rownames(original)
  candidate_rows <- rownames(candidate)
  if (!is.null(original_rows) && !is.null(candidate_rows)) {
    if (anyDuplicated(candidate_rows) || !setequal(original_rows, candidate_rows)) {
      stop(sprintf("%s returned incompatible sample row names.", method_name),
           call. = FALSE)
    }
    candidate <- candidate[original_rows, , drop = FALSE]
  }

  original_columns <- colnames(original)
  candidate_columns <- colnames(candidate)
  if (!is.null(original_columns) && !is.null(candidate_columns)) {
    if (anyDuplicated(candidate_columns) ||
        !setequal(original_columns, candidate_columns)) {
      stop(sprintf("%s returned incompatible protein column names.", method_name),
           call. = FALSE)
    }
    candidate <- candidate[, original_columns, drop = FALSE]
  }

  dimnames(candidate) <- dimnames(original)
  candidate
}

.trismart_resolve_workers <- function(
    workers = NULL,
    task_count = Inf,
    max_workers = 8L
) {
  max_workers <- .trismart_integer_scalar(
    max_workers,
    "max_workers",
    minimum = 1L
  )

  if (is.null(workers)) {
    detected <- parallel::detectCores(logical = TRUE)
    if (length(detected) != 1L || is.na(detected) || detected < 1L) {
      detected <- 1L
    }
    workers <- max(1L, detected - 1L)
  }

  workers <- .trismart_integer_scalar(workers, "workers", minimum = 1L)

  if (length(task_count) != 1L || is.na(task_count) || task_count < 0) {
    stop("task_count must be a non-negative number.", call. = FALSE)
  }
  if (is.finite(task_count)) {
    workers <- min(workers, max(1L, as.integer(task_count)))
  }
  min(workers, max_workers)
}

.trismart_copy_missing_only <- function(original, candidate, method_name) {
  candidate <- .trismart_align_candidate(original, candidate, method_name)

  out <- original
  missing <- is.na(original)
  out[missing] <- candidate[missing]

  if (anyNA(out)) {
    bad <- colnames(out)[colSums(is.na(out)) > 0L]
    stop(
      sprintf(
        "%s did not impute every requested missing value. Remaining columns: %s",
        method_name,
        paste(bad, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (any(!is.finite(out))) {
    stop(sprintf("%s produced non-finite values.", method_name), call. = FALSE)
  }

  dimnames(out) <- dimnames(original)
  out
}

# -----------------------------
# Zone classification
# -----------------------------

classify_missingness_zones <- function(
    X,
    protein_names = NULL,
    gmm_model = "VEE",
    residual_quantile = 0.90,
    seed = 123
) {
  .trismart_require("mgcv")
  .trismart_require("mclust")

  X <- .trismart_validate_matrix(X, min_rows = 3L, min_cols = 15L)

  if (!is.character(gmm_model) || length(gmm_model) != 1L ||
      is.na(gmm_model) || !nzchar(gmm_model)) {
    stop("gmm_model must be one non-empty Mclust model name.", call. = FALSE)
  }
  if (!is.numeric(residual_quantile) || length(residual_quantile) != 1L ||
      !is.finite(residual_quantile) || residual_quantile <= 0 ||
      residual_quantile >= 1) {
    stop("residual_quantile must be a single number strictly between 0 and 1.",
         call. = FALSE)
  }

  if (is.null(protein_names)) protein_names <- colnames(X)
  protein_names <- as.character(protein_names)
  if (length(protein_names) != ncol(X)) {
    stop("Length of protein_names must equal ncol(X).", call. = FALSE)
  }
  if (anyNA(protein_names) || any(trimws(protein_names) == "")) {
    stop("protein_names must be complete and non-empty.", call. = FALSE)
  }
  if (anyDuplicated(protein_names)) {
    stop("protein_names must be unique.", call. = FALSE)
  }

  median_abundance <- if (requireNamespace("matrixStats", quietly = TRUE)) {
    matrixStats::colMedians(X, na.rm = TRUE)
  } else {
    apply(X, 2L, stats::median, na.rm = TRUE)
  }
  missing_rate <- colSums(is.na(X)) / nrow(X)

  master <- data.frame(
    protein = protein_names,
    median_abund = as.numeric(median_abundance),
    miss_rate = as.numeric(missing_rate),
    usable = is.finite(median_abundance) & is.finite(missing_rate),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  usable_index <- which(master$usable)
  if (length(usable_index) < 15L) {
    stop("Too few valid proteins to model three missingness zones safely.",
         call. = FALSE)
  }
  model_data <- master[usable_index, , drop = FALSE]

  if (stats::sd(model_data$median_abund) <= sqrt(.Machine$double.eps)) {
    stop("Protein median abundance has zero variance; zone modeling is undefined.",
         call. = FALSE)
  }
  if (stats::sd(model_data$miss_rate) <= sqrt(.Machine$double.eps)) {
    stop("Protein missingness rate has zero variance; zone modeling is undefined.",
         call. = FALSE)
  }

  unique_abundance_count <- length(unique(model_data$median_abund))
  if (unique_abundance_count < 3L) {
    stop("At least three distinct protein median abundances are required.",
         call. = FALSE)
  }
  k_use <- min(
    10L,
    max(3L, floor(nrow(model_data) / 10L)),
    unique_abundance_count
  )

  fitted_objects <- .trismart_with_seed(seed, {
    s <- mgcv::s
    gam_fit <- mgcv::gam(
      miss_rate ~ s(median_abund, k = k_use),
      data = model_data,
      method = "REML"
    )

    prediction <- as.numeric(stats::predict(gam_fit, newdata = model_data))
    model_data$pred_miss_rate <- pmin(1, pmax(0, prediction))
    model_data$residual <- model_data$miss_rate - model_data$pred_miss_rate

    gmm_matrix <- scale(
      as.matrix(model_data[, c("median_abund", "miss_rate"), drop = FALSE])
    )
    if (any(!is.finite(gmm_matrix))) {
      stop("Scaling produced non-finite values for GMM fitting.", call. = FALSE)
    }

    gmm_fit <- tryCatch(
      mclust::Mclust(
        gmm_matrix,
        G = 3L,
        modelNames = gmm_model,
        verbose = FALSE
      ),
      error = function(e) {
        warning(
          sprintf(
            "GMM model '%s' failed (%s); using BIC-selected geometry with G = 3.",
            gmm_model,
            conditionMessage(e)
          ),
          call. = FALSE
        )
        mclust::Mclust(gmm_matrix, G = 3L, verbose = FALSE)
      }
    )

    list(
      gam_fit = gam_fit,
      gmm_fit = gmm_fit,
      model_data = model_data
    )
  })

  gam_fit <- fitted_objects$gam_fit
  gmm_fit <- fitted_objects$gmm_fit
  model_data <- fitted_objects$model_data

  if (is.null(gmm_fit) || is.null(gmm_fit$z) || is.null(gmm_fit$classification)) {
    stop("GMM clustering failed to return classifications and probabilities.",
         call. = FALSE)
  }

  posterior <- as.matrix(gmm_fit$z)
  clusters <- as.integer(gmm_fit$classification)
  cluster_ids <- sort(unique(clusters))
  if (length(cluster_ids) != 3L || ncol(posterior) != 3L) {
    stop("GMM did not produce exactly three non-empty clusters.", call. = FALSE)
  }
  model_data$gmm_cluster <- clusters

  cluster_summary <- do.call(
    rbind,
    lapply(cluster_ids, function(cluster_id) {
      rows <- clusters == cluster_id
      data.frame(
        cluster = cluster_id,
        median_abund = mean(model_data$median_abund[rows]),
        miss_rate = mean(model_data$miss_rate[rows]),
        stringsAsFactors = FALSE
      )
    })
  )

  safe_z <- function(x) {
    sx <- stats::sd(x)
    if (!is.finite(sx) || sx <= sqrt(.Machine$double.eps)) {
      rep(0, length(x))
    } else {
      (x - mean(x)) / sx
    }
  }

  cluster_summary$abund_z <- safe_z(cluster_summary$median_abund)
  cluster_summary$miss_z <- safe_z(cluster_summary$miss_rate)
  contrast <- cluster_summary$abund_z - cluster_summary$miss_z

  zone1_cluster <- cluster_summary$cluster[which.max(contrast)]
  zone3_cluster <- cluster_summary$cluster[which.min(contrast)]
  zone2_cluster <- setdiff(cluster_summary$cluster, c(zone1_cluster, zone3_cluster))

  if (zone1_cluster == zone3_cluster || length(zone2_cluster) != 1L) {
    warning(
      "Cluster contrast mapping was ambiguous; using abundance ordering.",
      call. = FALSE
    )
    abundance_order <- order(cluster_summary$median_abund, decreasing = TRUE)
    zone1_cluster <- cluster_summary$cluster[abundance_order[1L]]
    zone3_cluster <- cluster_summary$cluster[abundance_order[3L]]
    zone2_cluster <- cluster_summary$cluster[abundance_order[2L]]
  }

  posterior_column <- function(cluster_id) {
    if (!is.null(colnames(posterior)) &&
        as.character(cluster_id) %in% colnames(posterior)) {
      match(as.character(cluster_id), colnames(posterior))
    } else {
      match(cluster_id, cluster_ids)
    }
  }

  posterior_indices <- c(
    posterior_column(zone1_cluster),
    posterior_column(zone2_cluster),
    posterior_column(zone3_cluster)
  )
  if (anyNA(posterior_indices) || anyDuplicated(posterior_indices)) {
    stop("Could not map GMM posterior columns to the three zones.",
         call. = FALSE)
  }
  model_data$P_Zone1_GMM <- posterior[, posterior_indices[1L]]
  model_data$P_Zone2_GMM <- posterior[, posterior_indices[2L]]
  model_data$P_Zone3_GMM <- posterior[, posterior_indices[3L]]

  residual <- model_data$residual
  residual_threshold <- as.numeric(
    stats::quantile(abs(residual), residual_quantile, na.rm = TRUE, names = FALSE)
  )
  residual_iqr <- stats::IQR(residual, na.rm = TRUE)
  lambda <- if (!is.finite(residual_iqr) || residual_iqr <= 0) {
    1
  } else {
    log(2) / residual_iqr
  }

  delta_positive <- ifelse(
    residual > residual_threshold,
    1 - exp(-lambda * (residual - residual_threshold)),
    0
  )
  delta_negative <- ifelse(
    residual < -residual_threshold,
    1 - exp(-lambda * (abs(residual) - residual_threshold)),
    0
  )
  delta_positive <- pmin(1, pmax(0, delta_positive))
  delta_negative <- pmin(1, pmax(0, delta_negative))

  p1_final <- model_data$P_Zone1_GMM * (1 - delta_positive)
  p2_final <- model_data$P_Zone2_GMM +
    model_data$P_Zone1_GMM * delta_positive +
    model_data$P_Zone3_GMM * delta_negative
  p3_final <- model_data$P_Zone3_GMM * (1 - delta_negative)

  probability_sum <- p1_final + p2_final + p3_final
  if (any(!is.finite(probability_sum)) || any(probability_sum <= 0)) {
    stop("Zone probability refinement produced invalid probabilities.",
         call. = FALSE)
  }

  model_data$P_Zone1 <- p1_final / probability_sum
  model_data$P_Zone2 <- p2_final / probability_sum
  model_data$P_Zone3 <- p3_final / probability_sum

  probability_matrix <- as.matrix(
    model_data[, c("P_Zone1", "P_Zone2", "P_Zone3"), drop = FALSE]
  )
  model_data$final_zone <- c("Zone1", "Zone2", "Zone3")[
    max.col(probability_matrix, ties.method = "first")
  ]

  additional_columns <- c(
    "pred_miss_rate", "residual", "gmm_cluster",
    "P_Zone1_GMM", "P_Zone2_GMM", "P_Zone3_GMM",
    "P_Zone1", "P_Zone2", "P_Zone3", "final_zone"
  )
  final_table <- master
  for (column in additional_columns) {
    final_table[[column]] <- NA
    final_table[[column]][usable_index] <- model_data[[column]]
  }
  rownames(final_table) <- final_table$protein

  list(
    protein_zone_table = final_table,
    gam_fit = gam_fit,
    gmm_fit = gmm_fit,
    cluster_summary = cluster_summary,
    residual_threshold = residual_threshold,
    lambda = lambda
  )
}

# -----------------------------
# First-stage MAR imputation
# -----------------------------

impute_high_abundance <- function(
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
) {
  X_miss <- .trismart_validate_matrix(
    X_miss,
    min_rows = 2L,
    min_cols = 1L
  )

  impmethod <- match.arg(impmethod)
  knn_axis <- match.arg(knn_axis)

  rf_workers <- .trismart_integer_scalar(
    rf_workers,
    "rf_workers",
    minimum = 1L
  )
  rf_ntree <- .trismart_integer_scalar(rf_ntree, "rf_ntree", minimum = 1L)
  rf_maxiter <- .trismart_integer_scalar(
    rf_maxiter,
    "rf_maxiter",
    minimum = 1L
  )
  bpca_npcs <- .trismart_integer_scalar(
    bpca_npcs,
    "bpca_npcs",
    minimum = 1L
  )
  knn_k <- .trismart_integer_scalar(knn_k, "knn_k", minimum = 1L)
  knn_maxp <- .trismart_integer_scalar(
    knn_maxp,
    "knn_maxp",
    minimum = 1L
  )
  if (!is.numeric(knn_rowmax) || length(knn_rowmax) != 1L ||
      !is.finite(knn_rowmax) || knn_rowmax <= 0 || knn_rowmax > 1) {
    stop("knn_rowmax must be in (0, 1].", call. = FALSE)
  }
  if (!is.numeric(knn_colmax) || length(knn_colmax) != 1L ||
      !is.finite(knn_colmax) || knn_colmax <= 0 || knn_colmax > 1) {
    stop("knn_colmax must be in (0, 1].", call. = FALSE)
  }

  zone1_cols <- .trismart_validate_columns(zone1_cols, X_miss, "zone1_cols")
  zone2_cols <- .trismart_validate_columns(zone2_cols, X_miss, "zone2_cols")
  impute_cols <- unique(c(zone1_cols, zone2_cols))

  if (!length(impute_cols)) {
    return(X_miss)
  }

  X_sub <- X_miss[, impute_cols, drop = FALSE]
  if (!anyNA(X_sub)) return(X_sub)

  run_knn <- function() {
    .trismart_require("impute")
    .trismart_with_seed(seed, {
      knn_input <- X_sub 
      knn_result <- impute::impute.knn(
        knn_input,
        k = as.integer(knn_k),
        rowmax = knn_rowmax,
        colmax = knn_colmax,
        maxp = as.integer(knn_maxp),
        rng.seed = .trismart_seed_offset(seed)
      )$data
      })
  }

  run_bpca <- function() {
    .trismart_require("pcaMethods")
    .trismart_with_seed(seed, {
      n_pcs <- min(
        as.integer(bpca_npcs),
        nrow(X_sub) - 1L,
        ncol(X_sub) - 1L
      )
      if (n_pcs < 1L) {
        stop("BPCA requires at least two rows and two columns.", call. = FALSE)
      }
      result <- pcaMethods::pca(
        t(X_sub),
        method = "bpca",
        nPcs = n_pcs,
        center = TRUE,
        completeObs = TRUE
      )
      t(pcaMethods::completeObs(result))
    })
  }

  run_random_forest <- function() {
    .trismart_require("missForest")
    .trismart_with_seed(seed, {
      missforest_function <- getExportedValue("missForest", "missForest")
      missforest_formals <- names(formals(missforest_function))
      missforest_arguments <- list(
        xmis = as.data.frame(t(X_sub), check.names = FALSE),
        maxiter = as.integer(rf_maxiter),
        ntree = as.integer(rf_ntree),
        verbose = isTRUE(verbose)
      )

      # Newer ranger-backed missForest releases expose both arguments and
      # manage their own worker threads. Older releases with only
      # `parallelize` require an externally registered foreach backend, so
      # leave them at their safe sequential default.
     if ("num.threads" %in% missforest_formals) {
      missforest_arguments$num.threads <- 1L
      }

      if ("parallelize" %in% missforest_formals) {
      missforest_arguments$parallelize <- "no"
      }

      result <- do.call(missforest_function, missforest_arguments)
      as.matrix(t(result$ximp))
    })
  }

  candidate <- switch(
    impmethod,
    KNN = run_knn(),
    RandomForest = run_random_forest(),
    BPCA = run_bpca(),
    Ensemble = {
      knn_candidate <- run_knn()
      bpca_candidate <- run_bpca()
      rf_candidate <- run_random_forest()

      if (!identical(dim(knn_candidate), dim(X_sub)) ||
          !identical(dim(bpca_candidate), dim(X_sub)) ||
          !identical(dim(rf_candidate), dim(X_sub))) {
        stop("One or more ensemble members returned incompatible dimensions.",
             call. = FALSE)
      }
      (knn_candidate + bpca_candidate + rf_candidate) / 3
    }
  )

  .trismart_copy_missing_only(X_sub, candidate, impmethod)
}

# -----------------------------
# Second-stage mixture imputation
# -----------------------------

.trismart_select_predictors <- function(
    X_miss,
    X_imp,
    target_cols,
    predictor_pool,
    max_predictors,
    min_observed
) {
  predictor_list <- setNames(vector("list", length(target_cols)), target_cols)
  if (!length(target_cols) || !length(predictor_pool)) return(predictor_list)

  finite_predictors <- colSums(!is.finite(X_imp[, predictor_pool, drop = FALSE])) == 0L
  predictor_pool <- predictor_pool[finite_predictors]
  if (!length(predictor_pool)) return(predictor_list)

  observed_mask <- !is.na(X_miss[, target_cols, drop = FALSE])
  pattern_key <- apply(observed_mask, 2L, paste0, collapse = "")
  pattern_groups <- split(seq_along(target_cols), pattern_key)

  for (group_indices in pattern_groups) {
    observed_rows <- which(observed_mask[, group_indices[1L]])
    if (length(observed_rows) < min_observed) next

    group_targets <- target_cols[group_indices]
    x_group <- X_imp[observed_rows, predictor_pool, drop = FALSE]
    y_group <- X_miss[observed_rows, group_targets, drop = FALSE]

    correlation_matrix <- suppressWarnings(
      abs(stats::cor(x_group, y_group, use = "everything"))
    )
    correlation_matrix <- as.matrix(correlation_matrix)
    correlation_matrix[!is.finite(correlation_matrix)] <- 0

    for (local_index in seq_along(group_targets)) {
      target <- group_targets[local_index]
      correlations <- correlation_matrix[, local_index]

      self_index <- match(target, predictor_pool, nomatch = 0L)
      if (self_index > 0L) correlations[self_index] <- -Inf

      positive <- which(correlations > 0)
      if (!length(positive)) next

      keep_count <- min(max_predictors, length(positive))
      ordered <- positive[order(correlations[positive], decreasing = TRUE)]
      predictor_list[[target]] <- predictor_pool[ordered[seq_len(keep_count)]]
    }
  }

  predictor_list
}

.trismart_fit_one_mixture_target <- function(
    index,
    X_miss,
    X_imp,
    target_cols,
    predictor_list,
    fold_id_list,
    alpha,
    nlambda,
    min_observed
) {
  target <- target_cols[index]
  missing_rows <- which(is.na(X_miss[, target]))
  observed_rows <- which(!is.na(X_miss[, target]))

  if (!length(missing_rows) || length(observed_rows) < min_observed) {
    return(NULL)
  }

  predictors <- predictor_list[[target]]
  if (!length(predictors)) return(NULL)

  x_train <- X_imp[observed_rows, predictors, drop = FALSE]
  x_predict <- X_imp[missing_rows, predictors, drop = FALSE]
  y_train <- X_miss[observed_rows, target]

  if (any(!is.finite(x_train)) || any(!is.finite(x_predict)) ||
      any(!is.finite(y_train))) {
    return(NULL)
  }

  predictor_sd <- apply(x_train, 2L, stats::sd)
  usable_predictor <- is.finite(predictor_sd) & predictor_sd > sqrt(.Machine$double.eps)
  x_train <- x_train[, usable_predictor, drop = FALSE]
  x_predict <- x_predict[, usable_predictor, drop = FALSE]
  if (!ncol(x_train) || stats::sd(y_train) <= sqrt(.Machine$double.eps)) {
    return(NULL)
  }

  fold_id <- fold_id_list[[index]]
  fit <- tryCatch(
    glmnet::cv.glmnet(
      x = x_train,
      y = y_train,
      family = "gaussian",
      alpha = alpha,
      nlambda = as.integer(nlambda),
      foldid = fold_id,
      grouped = FALSE,
      standardize = TRUE,
      intercept = TRUE,
      parallel = FALSE
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)

  prediction <- tryCatch(
    as.numeric(stats::predict(fit, newx = x_predict, s = "lambda.min")),
    error = function(e) numeric(0)
  )
  if (length(prediction) != length(missing_rows) || any(!is.finite(prediction))) {
    return(NULL)
  }

  list(rows = missing_rows, values = prediction)
}

# This entry point is exported to PSOCK workers together with the objects below.
.trismart_mixture_worker_entry <- function(index) {
  .trismart_fit_one_mixture_target(
    index = index,
    X_miss = .trismart_worker_X_miss,
    X_imp = .trismart_worker_X_imp,
    target_cols = .trismart_worker_target_cols,
    predictor_list = .trismart_worker_predictor_list,
    fold_id_list = .trismart_worker_fold_id_list,
    alpha = .trismart_worker_alpha,
    nlambda = .trismart_worker_nlambda,
    min_observed = .trismart_worker_min_observed
  )
}

fast_mixture_impute <- function(
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
) {
  X_miss <- .trismart_validate_matrix(
    X_miss,
    min_rows = 3L,
    min_cols = 1L
  )
  X_imp_input <- .trismart_validate_matrix(
    X_impute_step1,
    min_rows = 3L,
    min_cols = 1L
  )

  mixture_idx <- .trismart_validate_columns(mixture_idx, X_miss, "mixture_idx")
  high_idx <- .trismart_validate_columns(high_idx, X_miss, "high_idx")
  required_stage1_columns <- unique(c(high_idx, mixture_idx))

  if (nrow(X_imp_input) != nrow(X_miss)) {
    stop("X_impute_step1 must have the same number of sample rows as X_miss.",
         call. = FALSE)
  }
  if (!identical(rownames(X_imp_input), rownames(X_miss))) {
    stop("X_impute_step1 and X_miss must have identical sample row names.",
         call. = FALSE)
  }
  unknown_stage1_columns <- setdiff(colnames(X_imp_input), colnames(X_miss))
  if (length(unknown_stage1_columns)) {
    stop(
      sprintf(
        "X_impute_step1 contains columns not found in X_miss: %s",
        paste(unknown_stage1_columns, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  missing_stage1_columns <- setdiff(
    required_stage1_columns,
    colnames(X_imp_input)
  )
  if (length(missing_stage1_columns)) {
    stop(
      sprintf(
        "X_impute_step1 is missing required high/mixture columns: %s",
        paste(missing_stage1_columns, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  # Accept both the original helper contract (a high/mixture subset) and a
  # full matrix. Rebuild one canonical full matrix by matching column names.
  X_imp <- X_miss
  if (length(required_stage1_columns)) {
    stage1_block <- X_imp_input[, required_stage1_columns, drop = FALSE]
    original_block <- X_miss[, required_stage1_columns, drop = FALSE]
    observed <- !is.na(original_block)
    if (!isTRUE(all.equal(
      stage1_block[observed],
      original_block[observed],
      tolerance = 0,
      check.attributes = FALSE
    ))) {
      stop("X_impute_step1 changed one or more originally observed values.",
           call. = FALSE)
    }
    if (anyNA(stage1_block) || any(!is.finite(stage1_block))) {
      stop("X_impute_step1 must completely and finitely impute high/mixture columns.",
           call. = FALSE)
    }
    X_imp[, required_stage1_columns] <- stage1_block
  }

  min_observed <- .trismart_integer_scalar(
    min_observed,
    "min_observed",
    minimum = 3L
  )
  nlambda <- .trismart_integer_scalar(nlambda, "nlambda", minimum = 2L)
  parallel_min_tasks <- .trismart_integer_scalar(
    parallel_min_tasks,
    "parallel_min_tasks",
    minimum = 1L
  )
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha < 0 || alpha > 1) {
    stop("alpha must be a single number between 0 and 1.", call. = FALSE)
  }

  missing_counts <- colSums(is.na(X_miss))
  target_cols <- intersect(mixture_idx, names(missing_counts)[missing_counts > 0L])

  X_final <- X_miss
  if (length(high_idx)) X_final[, high_idx] <- X_imp[, high_idx, drop = FALSE]

  if (!length(target_cols)) {
    result <- list(
      data = X_final,
      diagnostics = list(
        target_count = 0L,
        eligible_model_count = 0L,
        model_success_count = 0L,
        fallback_count = 0L,
        workers = 1L
      )
    )
    return(if (isTRUE(return_details)) result else result$data)
  }

  predictor_pool <- unique(c(high_idx, mixture_idx))
  if (is.null(max_predictors)) {
    max_predictors <- min(50L, max(10L, floor(nrow(X_miss) / 6L)))
  }
  max_predictors <- .trismart_integer_scalar(
    max_predictors,
    "max_predictors",
    minimum = 1L
  )

  predictor_list <- .trismart_select_predictors(
    X_miss = X_miss,
    X_imp = X_imp,
    target_cols = target_cols,
    predictor_pool = predictor_pool,
    max_predictors = max_predictors,
    min_observed = min_observed
  )

  observed_target_count <- colSums(!is.na(X_miss[, target_cols, drop = FALSE]))
  eligible_indices <- which(
    observed_target_count >= min_observed & lengths(predictor_list) > 0L
  )

  # If no target can support a model, use the documented Stage 1 fallback and
  # avoid loading glmnet or creating parallel workers.
  if (!length(eligible_indices)) {
    for (target in target_cols) {
      missing_rows <- which(is.na(X_miss[, target]))
      X_final[missing_rows, target] <- X_imp[missing_rows, target]
    }
    result <- list(
      data = X_final,
      diagnostics = list(
        target_count = length(target_cols),
        eligible_model_count = 0L,
        model_success_count = 0L,
        fallback_count = length(target_cols),
        workers = 1L,
        max_predictors = max_predictors
      )
    )
    return(if (isTRUE(return_details)) result else result$data)
  }

  .trismart_require("glmnet")

  fold_id_list <- vector("list", length(target_cols))
  generated_fold_ids <- .trismart_with_seed(seed, {
    lapply(eligible_indices, function(index) {
      n_observed <- observed_target_count[index]
      fold_count <- if (n_observed >= 35L) 5L else 3L
      fold_count <- min(fold_count, n_observed)
      sample(rep(seq_len(fold_count), length.out = n_observed))
    })
  })
  fold_id_list[eligible_indices] <- generated_fold_ids

  worker_count <- .trismart_resolve_workers(
    workers = workers,
    task_count = length(eligible_indices),
    max_workers = max_workers
  )
  use_parallel <- worker_count > 1L &&
    length(eligible_indices) >= parallel_min_tasks

  fit_sequentially <- function(indices) {
    lapply(indices, function(index) {
      .trismart_fit_one_mixture_target(
        index = index,
        X_miss = X_miss,
        X_imp = X_imp,
        target_cols = target_cols,
        predictor_list = predictor_list,
        fold_id_list = fold_id_list,
        alpha = alpha,
        nlambda = nlambda,
        min_observed = min_observed
      )
    })
  }

  eligible_results <- NULL
  if (use_parallel) {
    cluster <- NULL
    on.exit({
      if (!is.null(cluster)) {
        try(parallel::stopCluster(cluster), silent = TRUE)
      }
    }, add = TRUE)

    eligible_results <- tryCatch({
      cluster <- parallel::makeCluster(worker_count)

      .trismart_worker_X_miss <- X_miss
      .trismart_worker_X_imp <- X_imp
      .trismart_worker_target_cols <- target_cols
      .trismart_worker_predictor_list <- predictor_list
      .trismart_worker_fold_id_list <- fold_id_list
      .trismart_worker_alpha <- alpha
      .trismart_worker_nlambda <- nlambda
      .trismart_worker_min_observed <- min_observed

      parallel::clusterExport(
        cluster,
        varlist = c(
          ".trismart_fit_one_mixture_target",
          ".trismart_mixture_worker_entry",
          ".trismart_worker_X_miss",
          ".trismart_worker_X_imp",
          ".trismart_worker_target_cols",
          ".trismart_worker_predictor_list",
          ".trismart_worker_fold_id_list",
          ".trismart_worker_alpha",
          ".trismart_worker_nlambda",
          ".trismart_worker_min_observed"
        ),
        envir = environment()
      )
      parallel::clusterEvalQ(cluster, {
        if (!requireNamespace("glmnet", quietly = TRUE)) {
          stop("glmnet is not installed on a worker.")
        }
        NULL
      })

      result <- parallel::parLapplyLB(
        cluster,
        eligible_indices,
        .trismart_mixture_worker_entry
      )
      parallel::stopCluster(cluster)
      cluster <- NULL
      result
    }, error = function(e) {
      warning(
        sprintf(
          "Parallel mixture fitting failed (%s); retrying sequentially.",
          conditionMessage(e)
        ),
        call. = FALSE
      )
      if (!is.null(cluster)) {
        try(parallel::stopCluster(cluster), silent = TRUE)
        cluster <- NULL
      }
      NULL
    })
  }

  if (is.null(eligible_results)) {
    worker_count <- 1L
    eligible_results <- fit_sequentially(eligible_indices)
  }

  model_results <- setNames(vector("list", length(target_cols)), target_cols)
  model_results[eligible_indices] <- eligible_results
  success_count <- 0L
  for (target in target_cols) {
    result <- model_results[[target]]
    missing_rows <- which(is.na(X_miss[, target]))

    if (!is.null(result)) {
      X_final[result$rows, target] <- result$values
      success_count <- success_count + 1L
    } else {
      X_final[missing_rows, target] <- X_imp[missing_rows, target]
    }
  }

  if (anyNA(X_final[, unique(c(high_idx, mixture_idx)), drop = FALSE])) {
    bad <- unique(c(high_idx, mixture_idx))[
      colSums(is.na(X_final[, unique(c(high_idx, mixture_idx)), drop = FALSE])) > 0L
    ]
    stop(
      sprintf(
        "Stage 2 left missing values in high/mixture columns: %s",
        paste(bad, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  result <- list(
    data = X_final,
    diagnostics = list(
      target_count = length(target_cols),
      eligible_model_count = length(eligible_indices),
      model_success_count = success_count,
      fallback_count = length(target_cols) - success_count,
      workers = worker_count,
      max_predictors = max_predictors
    )
  )

  if (isTRUE(return_details)) result else result$data
}

# -----------------------------
# Third-stage MNAR imputation
# -----------------------------

# imputeLCMD follows the common omics convention: features/proteins in rows
# and samples in columns. TriSmart stores samples in rows and proteins in
# columns, so the default "proteins" axis applies the equivalent transposition.
# The "samples" option is retained only to reproduce the legacy script's direct
# calls, which treated sample rows as features.

.trismart_min_det_selected_proteins <- function(X, selected_columns, q = 0.01) {
  out <- X
  selected_index <- match(selected_columns, colnames(X))

  # Equivalent to t(imputeLCMD::impute.MinDet(t(X), q = q)), but only writes
  # the requested low-zone proteins and avoids creating another full matrix.
  low_quantile_by_sample <- apply(
    X,
    1L,
    stats::quantile,
    probs = q,
    na.rm = TRUE,
    names = FALSE
  )

  for (column_index in selected_index) {
    missing <- is.na(out[, column_index])
    out[missing, column_index] <- low_quantile_by_sample[missing]
  }
  out
}

.trismart_min_det_selected_legacy <- function(X, selected_columns, q = 0.01) {
  out <- X
  selected_index <- match(selected_columns, colnames(X))
  quantiles <- apply(
    X[, selected_index, drop = FALSE],
    2L,
    stats::quantile,
    probs = q,
    na.rm = TRUE,
    names = FALSE
  )

  for (j in seq_along(selected_index)) {
    column_index <- selected_index[j]
    missing <- is.na(out[, column_index])
    out[missing, column_index] <- quantiles[j]
  }
  out
}

.trismart_min_prob_selected_proteins <- function(
    X,
    selected_columns,
    q = 0.01,
    tune_sigma = 1
) {
  selected_index <- match(selected_columns, colnames(X))
  n_samples <- nrow(X)
  n_proteins <- ncol(X)

  # With proteins in rows and samples in columns, impute.MinProb estimates a
  # low quantile for each sample and a shared SD from protein-wise SDs.
  minimum_by_sample <- apply(
    X,
    1L,
    stats::quantile,
    probs = q,
    na.rm = TRUE,
    names = FALSE
  )

  sufficiently_observed <- colMeans(!is.na(X)) > 0.5
  if (!any(sufficiently_observed)) {
    stop("MinProb cannot estimate its shared standard deviation.", call. = FALSE)
  }
  protein_standard_deviation <- apply(
    X[, sufficiently_observed, drop = FALSE],
    2L,
    stats::sd,
    na.rm = TRUE
  )
  shared_sd <- stats::median(
    protein_standard_deviation,
    na.rm = TRUE
  ) * tune_sigma
  if (!is.finite(shared_sd) || shared_sd < 0) {
    stop("MinProb estimated an invalid standard deviation.", call. = FALSE)
  }

  # impute.MinProb(t(X)) draws n_proteins values for sample 1, then sample 2,
  # and so on. byrow = TRUE preserves that exact random-number ordering.
  draws <- matrix(
    stats::rnorm(
      n_samples * n_proteins,
      mean = rep(minimum_by_sample, each = n_proteins),
      sd = shared_sd
    ),
    nrow = n_samples,
    ncol = n_proteins,
    byrow = TRUE,
    dimnames = dimnames(X)
  )

  out <- X
  selected_block <- out[, selected_index, drop = FALSE]
  selected_draws <- draws[, selected_index, drop = FALSE]
  missing <- is.na(selected_block)
  selected_block[missing] <- selected_draws[missing]
  out[, selected_index] <- selected_block
  out
}

.trismart_min_prob_selected_legacy <- function(
    X,
    selected_columns,
    q = 0.01,
    tune_sigma = 1
) {
  selected_index <- match(selected_columns, colnames(X))
  n_samples <- nrow(X)

  minimum_by_column <- apply(
    X,
    2L,
    stats::quantile,
    probs = q,
    na.rm = TRUE,
    names = FALSE
  )

  sufficiently_observed <- rowMeans(!is.na(X)) > 0.5
  if (!any(sufficiently_observed)) {
    stop("MinProb cannot estimate its shared standard deviation.", call. = FALSE)
  }
  row_standard_deviation <- apply(
    X[sufficiently_observed, , drop = FALSE],
    1L,
    stats::sd,
    na.rm = TRUE
  )
  shared_sd <- stats::median(row_standard_deviation, na.rm = TRUE) * tune_sigma
  if (!is.finite(shared_sd) || shared_sd < 0) {
    stop("MinProb estimated an invalid standard deviation.", call. = FALSE)
  }

  draws <- matrix(
    stats::rnorm(
      length(X),
      mean = rep(minimum_by_column, each = n_samples),
      sd = shared_sd
    ),
    nrow = nrow(X),
    ncol = ncol(X),
    dimnames = dimnames(X)
  )

  out <- X
  selected_block <- out[, selected_index, drop = FALSE]
  selected_draws <- draws[, selected_index, drop = FALSE]
  missing <- is.na(selected_block)
  selected_block[missing] <- selected_draws[missing]
  out[, selected_index] <- selected_block
  out
}

.trismart_get_gsimp <- function() {
  if (exists("GS_impute", mode = "function", inherits = TRUE)) {
    return(get("GS_impute", mode = "function", inherits = TRUE))
  }
  if (requireNamespace("GSimp", quietly = TRUE) &&
      exists("GS_impute", envir = asNamespace("GSimp"), inherits = FALSE)) {
    return(get("GS_impute", envir = asNamespace("GSimp"), inherits = FALSE))
  }
  stop(
    paste0(
      "GSimp was requested, but GS_impute() is unavailable. ",
      "Install/load the GSimp implementation and all of its source files first."
    ),
    call. = FALSE
  )
}

impute_mnar_by_col <- function(
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
) {
  X_imp <- .trismart_validate_matrix(X_imp, min_rows = 2L, min_cols = 1L)
  method <- match.arg(method)
  imputelcmd_axis <- match.arg(imputelcmd_axis)
  if (identical(method, "Down-Shift")) {
    warning(
      "'Down-Shift' is retained as a backward-compatible alias for MinDet; it is not a Perseus-style Gaussian down-shift implementation.",
      call. = FALSE
    )
  }
  low_cols <- .trismart_validate_columns(low_cols, X_imp, "low_cols")

  if (!length(low_cols) || !anyNA(X_imp[, low_cols, drop = FALSE])) {
    return(X_imp)
  }
  if (!is.numeric(q) || length(q) != 1L || !is.finite(q) || q <= 0 || q >= 1) {
    stop("q must be strictly between 0 and 1.", call. = FALSE)
  }
  if (!is.numeric(tune_sigma) || length(tune_sigma) != 1L ||
      !is.finite(tune_sigma) || tune_sigma < 0) {
    stop("tune_sigma must be a finite non-negative number.", call. = FALSE)
  }

  gsimp_iters_each <- .trismart_integer_scalar(
    gsimp_iters_each,
    "gsimp_iters_each",
    minimum = 1L
  )
  gsimp_iters_all <- .trismart_integer_scalar(
    gsimp_iters_all,
    "gsimp_iters_all",
    minimum = 1L
  )
  gsimp_cores <- .trismart_integer_scalar(
    gsimp_cores,
    "gsimp_cores",
    minimum = 1L
  )

  original_missing <- is.na(X_imp)
  selected_index <- match(low_cols, colnames(X_imp))

  candidate <- .trismart_with_seed(seed, {
    switch(
      method,
      QRILC = {
        .trismart_require("imputeLCMD")
        if (imputelcmd_axis == "proteins") {
          completed <- imputeLCMD::impute.QRILC(
            t(X_imp),
            tune.sigma = tune_sigma
          )[[1L]]
          completed <- t(completed)
          dimnames(completed) <- dimnames(X_imp)
          completed
        } else {
          completed <- imputeLCMD::impute.QRILC(
            X_imp,
            tune.sigma = tune_sigma
          )[[1L]]
          dimnames(completed) <- dimnames(X_imp)
          completed
        }
      },
      MinProb = {
        if (imputelcmd_axis == "proteins") {
          .trismart_min_prob_selected_proteins(
            X_imp,
            selected_columns = low_cols,
            q = q,
            tune_sigma = tune_sigma
          )
        } else {
          .trismart_min_prob_selected_legacy(
            X_imp,
            selected_columns = low_cols,
            q = q,
            tune_sigma = tune_sigma
          )
        }
      },
      `Down-Shift` = {
        if (imputelcmd_axis == "proteins") {
          .trismart_min_det_selected_proteins(
            X_imp,
            selected_columns = low_cols,
            q = q
          )
        } else {
          .trismart_min_det_selected_legacy(
            X_imp,
            selected_columns = low_cols,
            q = q
          )
        }
      },
      MinDet = {
        if (imputelcmd_axis == "proteins") {
          .trismart_min_det_selected_proteins(
            X_imp,
            selected_columns = low_cols,
            q = q
          )
        } else {
          .trismart_min_det_selected_legacy(
            X_imp,
            selected_columns = low_cols,
            q = q
          )
        }
      },
      GSimp = {
        gsimp_function <- .trismart_get_gsimp()
        gsimp_result <- gsimp_function(
          data_miss = X_imp,
          iters_each = gsimp_iters_each,
          iters_all = gsimp_iters_all,
          initial = "qrilc",
          lo = -Inf,
          hi = "min",
          n_cores = gsimp_cores,
          imp_model = "glmnet_pred"
        )
        as.matrix(gsimp_result$data_imp)
      }
    )
  })

  candidate <- .trismart_align_candidate(X_imp, candidate, method)

  out <- X_imp
  selected_original <- X_imp[, selected_index, drop = FALSE]
  selected_candidate <- candidate[, selected_index, drop = FALSE]
  selected_missing <- is.na(selected_original)
  selected_original[selected_missing] <- selected_candidate[selected_missing]
  out[, selected_index] <- selected_original

  if (anyNA(out[, selected_index, drop = FALSE])) {
    bad <- low_cols[colSums(is.na(out[, selected_index, drop = FALSE])) > 0L]
    stop(
      sprintf(
        "%s left missing values in low-abundance columns: %s",
        method,
        paste(bad, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (any(!is.finite(out[, selected_index, drop = FALSE]))) {
    stop(sprintf("%s produced non-finite imputations.", method), call. = FALSE)
  }

  outside_selected <- setdiff(seq_len(ncol(X_imp)), selected_index)
  if (length(outside_selected) &&
      !identical(out[, outside_selected, drop = FALSE],
                 X_imp[, outside_selected, drop = FALSE])) {
    stop("Internal error: MNAR imputation modified a non-low-zone column.",
         call. = FALSE)
  }

  observed_selected <- !original_missing[, selected_index, drop = FALSE]
  if (!isTRUE(all.equal(
    out[, selected_index, drop = FALSE][observed_selected],
    X_imp[, selected_index, drop = FALSE][observed_selected],
    tolerance = 0,
    check.attributes = FALSE
  ))) {
    stop("Internal error: MNAR imputation changed observed values.",
         call. = FALSE)
  }

  dimnames(out) <- dimnames(X_imp)
  out
}

# -----------------------------
# End-to-end TriSmart pipeline
# -----------------------------

TriSmart <- function(
    X_miss,
    firstImpute,
    secondImpute,
    workers = NULL,
    max_workers = 8L,
    seed = 1L,
    residual_quantile = 0.90,
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
) {
  X_miss <- .trismart_validate_matrix(
    X_miss,
    min_rows = 3L,
    min_cols = 15L
  )

  firstImpute <- match.arg(
    as.character(firstImpute),
    c("KNN", "RandomForest", "BPCA", "Ensemble")
  )
  secondImpute <- match.arg(
    as.character(secondImpute),
    c("QRILC", "MinProb", "Down-Shift", "MinDet", "GSimp")
  )
  knn_axis <- match.arg(knn_axis)
  imputelcmd_axis <- match.arg(imputelcmd_axis)

  resolved_workers <- .trismart_resolve_workers(
    workers = workers,
    max_workers = max_workers
  )
  gsimp_cores <- .trismart_integer_scalar(
    gsimp_cores,
    "gsimp_cores",
    minimum = 1L
  )

  if (isTRUE(verbose)) {
    message(
      sprintf(
        "TriSmart: %d samples, %d proteins, %.2f%% missing, up to %d workers.",
        nrow(X_miss),
        ncol(X_miss),
        100 * mean(is.na(X_miss)),
        resolved_workers
      )
    )
  }

  zone_result <- classify_missingness_zones(
    X = X_miss,
    gmm_model = gmm_model,
    residual_quantile = residual_quantile,
    seed = .trismart_seed_offset(seed, 0L)
  )
  zone_table <- zone_result$protein_zone_table

  if (anyNA(zone_table$final_zone)) {
    bad <- zone_table$protein[is.na(zone_table$final_zone)]
    stop(
      sprintf("Some proteins could not be assigned to a zone: %s",
              paste(bad, collapse = ", ")),
      call. = FALSE
    )
  }

  high_idx <- zone_table$protein[zone_table$final_zone == "Zone1"]
  mixture_idx <- zone_table$protein[zone_table$final_zone == "Zone2"]
  low_idx <- zone_table$protein[zone_table$final_zone == "Zone3"]

  zone_union <- c(high_idx, mixture_idx, low_idx)
  if (length(zone_union) != ncol(X_miss) || anyDuplicated(zone_union) ||
      !setequal(zone_union, colnames(X_miss))) {
    stop("Zone assignment did not form a complete, disjoint partition.",
         call. = FALSE)
  }

  if (isTRUE(verbose)) {
    message(
      sprintf(
        "Zones: Zone1=%d, Zone2=%d, Zone3=%d.",
        length(high_idx), length(mixture_idx), length(low_idx)
      )
    )
  }

  high_mixture_columns <- unique(c(high_idx, mixture_idx))
  stage1_subset <- impute_high_abundance(
    X_miss = X_miss,
    zone1_cols = high_idx,
    zone2_cols = mixture_idx,
    impmethod = firstImpute,
    rf_ntree = rf_ntree,
    rf_maxiter = rf_maxiter,
    bpca_npcs = bpca_npcs,
    knn_axis = knn_axis,
    rf_workers = resolved_workers,
    seed = .trismart_seed_offset(seed, 1L),
    verbose = verbose
  )

  stage1 <- X_miss
  if (length(high_mixture_columns)) {
    stage1[, high_mixture_columns] <- stage1_subset
  }

  if (anyNA(stage1[, high_mixture_columns, drop = FALSE])) {
    bad <- high_mixture_columns[
      colSums(is.na(stage1[, high_mixture_columns, drop = FALSE])) > 0L
    ]
    stop(
      sprintf("Stage 1 left missing values in: %s", paste(bad, collapse = ", ")),
      call. = FALSE
    )
  }

  stage2_result <- fast_mixture_impute(
    X_miss = X_miss,
    X_impute_step1 = stage1,
    mixture_idx = mixture_idx,
    high_idx = high_idx,
    min_observed = min_observed,
    max_predictors = max_predictors,
    alpha = glmnet_alpha,
    nlambda = glmnet_nlambda,
    workers = resolved_workers,
    max_workers = max_workers,
    seed = .trismart_seed_offset(seed, 2L),
    return_details = TRUE
  )
  stage2 <- stage2_result$data

  stage3 <- impute_mnar_by_col(
    X_imp = stage2,
    low_cols = low_idx,
    method = secondImpute,
    q = mnar_q,
    tune_sigma = mnar_tune_sigma,
    gsimp_cores = gsimp_cores,
    imputelcmd_axis = imputelcmd_axis,
    seed = .trismart_seed_offset(seed, 3L)
  )

  if (anyNA(stage3)) {
    bad <- colnames(stage3)[colSums(is.na(stage3)) > 0L]
    stop(
      sprintf("Final output still contains missing values in: %s",
              paste(bad, collapse = ", ")),
      call. = FALSE
    )
  }
  if (any(!is.finite(stage3))) {
    stop("Final output contains non-finite values.", call. = FALSE)
  }

  observed <- !is.na(X_miss)
  if (!isTRUE(all.equal(
    stage3[observed],
    X_miss[observed],
    tolerance = 0,
    check.attributes = FALSE
  ))) {
    stop("Internal error: one or more observed values were changed.",
         call. = FALSE)
  }

  dimnames(stage3) <- dimnames(X_miss)

  details <- list(
    imputed_data = stage3,
    zone_result = zone_result,
    zone_columns = list(
      Zone1 = high_idx,
      Zone2 = mixture_idx,
      Zone3 = low_idx
    ),
    mixture_diagnostics = stage2_result$diagnostics,
    settings = list(
      firstImpute = firstImpute,
      secondImpute = secondImpute,
      workers = resolved_workers,
      seed = seed,
      knn_axis = knn_axis,
      imputelcmd_axis = imputelcmd_axis,
      rf_workers = resolved_workers,
      gsimp_cores = gsimp_cores
    )
  )

  if (isTRUE(return_details)) details else details$imputed_data
}

# Lower-case alias for projects that follow snake_case naming conventions.
tri_smart <- TriSmart
