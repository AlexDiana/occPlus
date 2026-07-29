
logistic <- function(x){
  1 / (1 + exp(-x))
}

reparamFactorModel <- function(U_output, L_output){

  L_output_reparam <- L_output
  U_output_reparam <- U_output

  d <- dim(L_output)[1]
  nchain <- dim(L_output)[4]
  niter <- dim(L_output)[3]

  for (chain in 1:nchain) {

    for (iter in 1:niter) {

      if(d == 1){

        L1 <- L_output[1,1,iter,chain]
        L_output_reparam[1,,iter,chain] <- L_output[1,,iter,chain] / L1
        U_output_reparam[,1,iter,chain] <- U_output[,1,iter,chain] * L1

      } else {

        L_current <- L_output[,,iter,chain]
        U_current <- U_output[,,iter,chain]

        qr_decomp <- qr(L_current)
        Q_current <- qr.Q(qr_decomp)
        R_current <- qr.R(qr_decomp)

        Q2 <- Q_current %*% diag(diag(R_current), nrow = d)
        invQ2 <- diag(1 / diag(R_current), nrow = d) %*% t(Q_current)

        L_new <- invQ2 %*% L_current
        U_new <- U_current %*% Q2

        L_output_reparam[,,iter,chain] <- L_new
        U_output_reparam[,,iter,chain] <- U_new
      }

    }
  }

  list("L_output" = L_output_reparam,
       "U_output" = U_output_reparam)
}

# COVARIATES AND SPLINES FUNCTIONS --------

create_covariates_matrix <- function(df,
                                     spline_vars = TRUE, # true for all or the names
                                     df_spline = 5,
                                     remove_intercept = TRUE,
                                     min_obs_per_df = 10) {

  if (is.matrix(df)) df <- as.data.frame(df)

  if (any(is.na(df))) stop("NA in covariates matrix")

  names_df   <- colnames(df)
  df         <- df %>% dplyr::mutate(dplyr::across(dplyr::where(~ !is.numeric(.x)), as.factor))
  is_numeric <- sapply(df, is.numeric)

  # Determine which specific numeric variables get splines
  if (is.logical(spline_vars)) {
    if (spline_vars) {
      target_spline_vars <- names_df[is_numeric] # All numeric
    } else {
      target_spline_vars <- character(0)          # None
    }
  } else if (is.character(spline_vars)) {
    # Check that requested variables actually exist and are numeric
    invalid_vars <- setdiff(spline_vars, names_df)
    if (length(invalid_vars) > 0) {
      stop(paste("Requested spline_vars not found in data:", paste(invalid_vars, collapse = ", ")))
    }
    non_num <- spline_vars[!is_numeric[spline_vars]]
    if (length(non_num) > 0) {
      warning(paste("Ignoring non-numeric variables in spline_vars:", paste(non_num, collapse = ", ")))
    }
    target_spline_vars <- intersect(spline_vars, names_df[is_numeric])
  } else {
    stop("'spline_vars' must be a logical (TRUE/FALSE) or a character vector of column names.")
  }

  mean_df    <- rep(NA, ncol(df))
  sd_df      <- rep(NA, ncol(df))
  cat_levels <- vector("list", ncol(df))
  bs_info    <- vector("list", ncol(df))  # Stores spline attributes

  names(mean_df) <- names(sd_df) <- names(cat_levels) <- names(bs_info) <- names_df

  # Standardize numericals and record factor levels
  for (col in 1:ncol(df)) {
    col_name <- names_df[col]
    if (is_numeric[col]) {
      mean_df[col] <- mean(df[[col]], na.rm = TRUE)
      sd_df[col]   <- sd(df[[col]], na.rm = TRUE)
      df[[col]]    <- (df[[col]] - mean_df[col]) / sd_df[col]
    } else {
      cat_levels[[col]] <- levels(df[[col]])
    }
  }

  # Sample size check rule for splines
  min_required_obs <- df_spline * min_obs_per_df

  # Build matrix columns
  list_X <- lapply(1:ncol(df), function(col) {

    col_name <- names_df[col]
    col_data <- df[[col]]

    # Check if variable is selected for spline AND meets sample size rule
    should_spline <- (col_name %in% target_spline_vars) &&
      (length(unique(col_data)) >= min_required_obs)

    if (col_name %in% target_spline_vars && length(unique(col_data)) < min_required_obs) {
      message(paste0("Note: '", col_name, "' has fewer than ", min_required_obs,
                     " unique values. Fitting as linear term instead of spline."))
    }

    if (is_numeric[col] && should_spline) {

      # Generate spline basis
      bs_obj <- splines::bs(col_data, df = df_spline, intercept = FALSE)

      # Save knot info for prediction step later
      bs_info[[col]] <<- list(
        knots          = attr(bs_obj, "knots"),
        Boundary.knots = attr(bs_obj, "Boundary.knots"),
        degree         = attr(bs_obj, "degree"),
        intercept      = attr(bs_obj, "intercept")
      )

      colnames(bs_obj) <- paste0(col_name, "_s", seq_len(ncol(bs_obj)))
      return(bs_obj)

    } else if (is_numeric[col]) {
      # Linear term fallback
      mat <- matrix(col_data, ncol = 1)
      colnames(mat) <- col_name
      return(mat)

    } else {
      # Categorical variables
      return(data.frame(col_data, fix.empty.names = FALSE) %>% setNames(col_name))
    }
  })

  df_transformed <- do.call(cbind, list_X)

  # Generate dummy design matrix
  X_final <- stats::model.matrix(~ ., data = as.data.frame(df_transformed))

  if (remove_intercept) {
    X_final <- X_final[, -1, drop = FALSE]
  }

  list_matrix <- list(
    names_df           = names_df,
    mean_df            = mean_df,
    sd_df              = sd_df,
    cat_levels         = cat_levels,
    is_numeric         = is_numeric,
    bs_info            = bs_info,
    target_spline_vars = target_spline_vars
  )

  return(list(X = X_final, list_matrix = list_matrix, X0 = df))
}

transform_new_covariates <- function(df, list_matrix, remove_intercept = TRUE) {

  names_df   <- list_matrix$names_df
  mean_df    <- list_matrix$mean_df
  sd_df      <- list_matrix$sd_df
  cat_levels <- list_matrix$cat_levels
  is_numeric <- list_matrix$is_numeric
  bs_info    <- list_matrix$bs_info

  if (is.matrix(df)) df <- as.data.frame(df)

  if (any(is.na(df))) stop("NA in covariates matrix")

  missing_cols <- setdiff(names_df, colnames(df))
  if (length(missing_cols) > 0) {
    stop(paste("New data is missing required columns:", paste(missing_cols, collapse = ", ")))
  }

  df <- df[, names_df, drop = FALSE] %>%
    dplyr::mutate(dplyr::across(dplyr::where(~ !is.numeric(.x)), as.factor))

  # 1. Standardize numericals and enforce factor levels
  for (col in 1:ncol(df)) {
    if (is_numeric[col]) {
      df[, col] <- (df[, col] - mean_df[col]) / sd_df[col]
    } else {
      df[[col]] <- factor(df[[col]], levels = cat_levels[[col]])
    }
  }

  # 2. Evaluate splines or keep linear
  list_X <- lapply(1:ncol(df), function(col) {
    col_name <- names_df[col]
    col_data <- df[[col]]

    if (is_numeric[col] && !is.null(bs_info[[col]])) {
      # Reconstruct spline using saved training knots
      info <- bs_info[[col]]
      bs_obj <- splines::bs(
        col_data,
        knots          = info$knots,
        Boundary.knots = info$Boundary.knots,
        degree         = info$degree,
        intercept      = info$intercept
      )

      colnames(bs_obj) <- paste0(col_name, "_s", seq_len(ncol(bs_obj)))
      return(bs_obj)

    } else if (is_numeric[col]) {
      mat <- matrix(col_data, ncol = 1)
      colnames(mat) <- col_name
      return(mat)

    } else {
      return(data.frame(col_data, fix.empty.names = FALSE) %>% setNames(col_name))
    }
  })

  df_transformed <- do.call(cbind, list_X)

  X_new <- stats::model.matrix(~ ., data = as.data.frame(df_transformed))

  if (remove_intercept) {
    X_new <- X_new[, -1, drop = FALSE]
  }

  return(X_new)
}

transformCovariatesMatrix <- function(df, list_matrix, remove_intercept){

  names_df <- list_matrix$names_df
  mean_df <- list_matrix$mean_df
  sd_df <- list_matrix$sd_df
  cat_levels <- list_matrix$cat_levels
  is_numeric <- list_matrix$is_numeric

  if(any(is.na(df))) stop("NA in covariates matrix")

  # missing columns
  missing_cols <- setdiff(names_df, colnames(df))
  if (length(missing_cols) > 0) {
    stop(paste("New data is missing critical columns required by the model:",
               paste(missing_cols, collapse = ", ")))
  }

  # sort the new data and convert non numeric to
  df <- df[,names_df,drop=F] %>%
    dplyr::mutate(dplyr::across(dplyr::where(~ !is.numeric(.x)), as.factor))

  # standardise numerical and categorical
  for (col in 1:ncol(df)) {
    if(is_numeric[col]){

      df[,col] <- (df[,col] - mean_df[col]) / sd_df[col]

    }else{

      levels(df[[col]]) <- cat_levels[[col]]

    }
  }

  df <- stats::model.matrix(~., df)

  # Remove intercept if requested (first column is always intercept)
  if (remove_intercept) {
    df <- df[, -1, drop = FALSE]
  }

  return(df)

}

createSplinesObjects <- function(X, df){

  list_ns <- lapply(seq_len(ncol(X)), function(j) {
    ns_j <- bs(X[, j], df = 5, intercept = F)
    ns_j
  })

  list_ns

}

createSplinesMatrixSingleCov <- function(n_s, X_cov){

  Zj <- predict(n_s, X_cov)
  colnames(Zj) <- paste0(colnames(X_cov), " - s", seq_len(ncol(Zj)))

  Zj
}

createSplinesMatrix <- function(list_ns, X_new){

  Z_list <- lapply(seq_len(ncol(X_new)), function(j) {
    # Zj <- predict(list_ns[[j]], X_new[,j])
    # colnames(Zj) <- paste0(colnames(X)[j], " - s", seq_len(ncol(Zj)))
    Zj <- createSplinesMatrixSingleCov(list_ns[[j]], X_new[,j,drop=F])
    Zj
  })

  Z <- do.call(cbind, Z_list)

  Z

}

returnCovariateEffect_base <- function(cov_name,
                                       idx_species,
                                       sp_name,
                                       B0_output_vec,
                                       B_output_vec,
                                       list_matrix,
                                       speciesNames,
                                       X0, X,
                                       n_points = 200,
                                       link = c("identity", "logit"),
                                       confidence = .95){

  conflevels <- c((1 - confidence)/2, .5, (1 + confidence)/2)

  # prepare covariates info
  {
    is_num    <- list_matrix$is_numeric[[cov_name]]
    is_spline <- !is.null(list_matrix$bs_info[[cov_name]])

    if (is_num) {
      # Sequence on raw scale
      cov_min <- min(X0[[cov_name]], na.rm = TRUE)
      cov_max <- max(X0[[cov_name]], na.rm = TRUE)
      cov_seq_raw <- seq(cov_min, cov_max, length.out = n_points)
      cov_std <- (cov_seq_raw - list_matrix$mean_df[[cov_name]]) / list_matrix$sd_df[[cov_name]]

      if (is_spline) {
        # Re-apply splines using stored knots on standardized data
        info <- list_matrix$bs_info[[cov_name]]
        X_sub <- splines::bs(
          cov_std,
          knots          = info$knots,
          Boundary.knots = info$Boundary.knots,
          degree         = info$degree,
          intercept      = info$intercept
        )
        colnames(X_sub) <- paste0(cov_name, "_s", seq_len(ncol(X_sub)))

      } else {
        # Linear numeric term
        X_sub <- matrix(cov_std, ncol = 1)
        colnames(X_sub) <- cov_name
      }

    } else {

      all_levels     <- list_matrix$cat_levels[[cov_name]]
      baseline_level <- all_levels[1]
      active_levels  <- all_levels[-1]

      temp_df <- data.frame(val = factor(active_levels, levels = all_levels))
      colnames(temp_df) <- cov_name

      temp_X <- stats::model.matrix(~ ., data = temp_df)
      X_sub  <- temp_X[, -1, drop = FALSE]
    }

    cov_indices <- sapply(colnames(X_sub), function(name){
      grep(name, colnames(X))
    })

    if (any(is.na(cov_indices))) {
      stop(paste("Could not find exact columns in MCMC output for:", cov_name))
    }

  }

  all_species_data <- data.frame()

  for (i in seq_along(idx_species)) {

    sp_idx  <- idx_species[i]
    # Index by sp_idx, not the loop counter i -- speciesNames[i] mislabels
    # every species whenever idx_species isn't the prefix 1:k (same defect
    # as TODO.md group C item 1, found here as part of item 5).
    sp_name <- speciesNames[sp_idx]

    beta_mcmc_j <- B_output_vec[, , sp_idx]
    beta_mcmc_sub <- beta_mcmc_j[, cov_indices, drop = FALSE]

    partial_effect <- X_sub %*% t(beta_mcmc_sub)

    if(link == "logit") partial_effect <- logistic(partial_effect)

    # 6. Format data for ggplot depending on variable type
    if (is_num) {

      intercept_draws <- B0_output_vec[, sp_idx]
      total_effect    <- sweep(partial_effect, 2, intercept_draws, FUN = "+")

      # Numeric: Summarize to mean and 95% Credible Intervals
      sp_data <- data.frame(
        x       = cov_seq_raw,
        mean    = apply(total_effect, 1, quantile, probs = conflevels[2]),
        lower   = apply(total_effect, 1, quantile, probs = conflevels[1]),
        upper   = apply(total_effect, 1, quantile, probs = conflevels[3]),
        Species = sp_name
      )

    } else {
      # Categorical: Keep all draws to generate boxplots
      # Transpose so rows are grid points, columns are draws, then pivot long
      sp_data <- as.data.frame(partial_effect) %>%
        mutate(x = factor(active_levels, levels = active_levels), Species = sp_name) %>%
        pivot_longer(cols = -c(x, Species), names_to = "draw", values_to = "value")
    }


    all_species_data <- bind_rows(all_species_data, sp_data)

  }

  all_species_data

}

plotCovariateEffect_base <- function(idx_species,
                               cov_names,
                               B0_output_vec,
                               B_output_vec,
                               list_matrix,
                               speciesNames,
                               X0, X,
                               n_points = 200,
                               link = c("identity", "logit"),
                               confidence = .95) {

  X0 <- as.data.frame(X0)
  X <- as.data.frame(X)

  # B0_output_vec/B_output_vec arrive already collapsed to (draws x species)
  # and (draws x covariate x species) by the caller (plotCovariateEffect()).
  # A previous version re-applied apply(..., c(1,2), c) here, which collapsed
  # the species margin a second time and left B_output_vec's third dimension
  # sized by ncov_psi instead of S -- "subscript out of bounds" for any
  # sp_idx > ncov_psi. Fixed as part of TODO.md group C item 5.

  plot_list <- list()

  for (cov_name in cov_names) {

    is_num    <- list_matrix$is_numeric[[cov_name]]
    is_spline <- !is.null(list_matrix$bs_info[[cov_name]])

    # prepare covariates info
    {
      is_num    <- list_matrix$is_numeric[[cov_name]]
      is_spline <- !is.null(list_matrix$bs_info[[cov_name]])

      if (is_num) {
        # Sequence on raw scale
        cov_min <- min(X0[[cov_name]], na.rm = TRUE)
        cov_max <- max(X0[[cov_name]], na.rm = TRUE)
        cov_seq_raw <- seq(cov_min, cov_max, length.out = n_points)
        cov_std <- (cov_seq_raw - list_matrix$mean_df[[cov_name]]) / list_matrix$sd_df[[cov_name]]

        if (is_spline) {
          # Re-apply splines using stored knots on standardized data
          info <- list_matrix$bs_info[[cov_name]]
          X_sub <- splines::bs(
            cov_std,
            knots          = info$knots,
            Boundary.knots = info$Boundary.knots,
            degree         = info$degree,
            intercept      = info$intercept
          )
          colnames(X_sub) <- paste0(cov_name, "_s", seq_len(ncol(X_sub)))

        } else {
          # Linear numeric term
          X_sub <- matrix(cov_std, ncol = 1)
          colnames(X_sub) <- cov_name
        }

      } else {

        all_levels     <- list_matrix$cat_levels[[cov_name]]
        baseline_level <- all_levels[1]
        active_levels  <- all_levels[-1]

        temp_df <- data.frame(val = factor(active_levels, levels = all_levels))
        colnames(temp_df) <- cov_name

        temp_X <- stats::model.matrix(~ ., data = temp_df)
        X_sub  <- temp_X[, -1, drop = FALSE]
      }

      cov_indices <- sapply(colnames(X_sub), function(name){
        grep(name, colnames(X))
      })

      if (any(is.na(cov_indices))) {
        stop(paste("Could not find exact columns in MCMC output for:", cov_name))
      }

    }

    if(F){
      all_species_data <- data.frame()

      for (i in seq_along(idx_species)) {

        sp_idx  <- idx_species[i]
        sp_name <- speciesNames[i]

        # prepare covariates info
        if(F){
          is_num    <- list_matrix$is_numeric[[cov_name]]
          is_spline <- !is.null(list_matrix$bs_info[[cov_name]])

          if (is_num) {
            # Sequence on raw scale
            cov_min <- min(X0[[cov_name]], na.rm = TRUE)
            cov_max <- max(X0[[cov_name]], na.rm = TRUE)
            cov_seq_raw <- seq(cov_min, cov_max, length.out = n_points)
            cov_std <- (cov_seq_raw - list_matrix$mean_df[[cov_name]]) / list_matrix$sd_df[[cov_name]]

            if (is_spline) {
              # Re-apply splines using stored knots on standardized data
              info <- list_matrix$bs_info[[cov_name]]
              X_sub <- splines::bs(
                cov_std,
                knots          = info$knots,
                Boundary.knots = info$Boundary.knots,
                degree         = info$degree,
                intercept      = info$intercept
              )
              colnames(X_sub) <- paste0(cov_name, "_s", seq_len(ncol(X_sub)))

            } else {
              # Linear numeric term
              X_sub <- matrix(cov_std, ncol = 1)
              colnames(X_sub) <- cov_name
            }

          } else {

            all_levels     <- list_matrix$cat_levels[[cov_name]]
            baseline_level <- all_levels[1]
            active_levels  <- all_levels[-1]

            temp_df <- data.frame(val = factor(active_levels, levels = all_levels))
            colnames(temp_df) <- cov_name

            temp_X <- stats::model.matrix(~ ., data = temp_df)
            X_sub  <- temp_X[, -1, drop = FALSE]
          }

          cov_indices <- sapply(colnames(X_sub), function(name){
            grep(name, colnames(X))
          })

          if (any(is.na(cov_indices))) {
            stop(paste("Could not find exact columns in MCMC output for:", cov_name))
          }

        }

        # species info
        if(F){

          beta_mcmc_j <- B_output_vec[, , sp_idx]
          beta_mcmc_sub <- beta_mcmc_j[, cov_indices, drop = FALSE]

          partial_effect <- X_sub %*% t(beta_mcmc_sub)

          if(link == "logit") partial_effect <- logistic(partial_effect)

          # 6. Format data for ggplot depending on variable type
          if (is_num) {

            intercept_draws <- B0_output_vec[, sp_idx]
            total_effect    <- sweep(partial_effect, 2, intercept_draws, FUN = "+")

            # Numeric: Summarize to mean and 95% Credible Intervals
            sp_data <- data.frame(
              x       = cov_seq_raw,
              mean    = apply(total_effect, 1, mean),
              lower   = apply(total_effect, 1, quantile, probs = 0.025),
              upper   = apply(total_effect, 1, quantile, probs = 0.975),
              Species = sp_name
            )
          } else {
            # Categorical: Keep all draws to generate boxplots
            # Transpose so rows are grid points, columns are draws, then pivot long
            sp_data <- as.data.frame(partial_effect) %>%
              mutate(x = factor(active_levels, levels = active_levels), Species = sp_name) %>%
              pivot_longer(cols = -c(x, Species), names_to = "draw", values_to = "value")
          }
        }

        sp_data <- returnCovariateEffect_base(
          cov_name,
          sp_idx,
          sp_name,
          B0_output_vec,
          B_output_vec,
          list_matrix,
          speciesNames,
          X0, X,
          n_points = n_points,
          link = link,
          confidence
        )

        all_species_data <- bind_rows(all_species_data, sp_data)
      }

    }

    all_species_data <- returnCovariateEffect_base(
      cov_name,
      idx_species,
      sp_name,
      B0_output_vec,
      B_output_vec,
      list_matrix,
      speciesNames,
      X0, X,
      n_points,
      link = link,
      confidence
    )

    # 7. Generate the Plot
    if (is_num) {
      type_title <- if (is_spline) "Spline" else "Linear"

      # Automatically pick ~5 clean tick values across the raw range
      raw_ticks <- pretty(cov_seq_raw, n = 5)

      # Plot for continuous (Line + Ribbon)
      p <- ggplot(all_species_data, aes(x = x, y = mean)) +
        geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#3388ff", alpha = 0.3) +
        geom_line(color = "#0044cc", linewidth = 1) +
        scale_x_continuous(breaks = raw_ticks) +
        facet_wrap(~ Species, scales = "free_y") +
        labs(title = paste("Effect of", cov_name), x = cov_name, y = "Linear Predictor") +
        theme_bw()
    } else {
      # Plot for categorical (Boxplot across MCMC draws)
      p <- ggplot(all_species_data, aes(x = factor(x), y = value)) +
        geom_boxplot(fill = "#e6f2ff", color = "#0044cc", outlier.alpha = 0.1) +
        facet_wrap(~ Species, scales = "free_y") +
        labs(title = paste("Effect of", cov_name), x = cov_name, y = "Linear Predictor") +
        theme_bw()
    }

    # Store plot in list using covariate name
    plot_list[[cov_name]] <- p
  }

  return(plot_list)
}

# SPATIAL FUNCTIONS -----------

buildGrid <- function(XY_sp, gridStep){

  x_grid <- seq(min(XY_sp[,1]) - (1.5) * gridStep,
                max(XY_sp[,1]) + (1.5) * gridStep, by = gridStep)
  y_grid <- seq(min(XY_sp[,2]) - (1.5) * gridStep,
                max(XY_sp[,2]) + (1.5) * gridStep, by = gridStep)

  pointInGrid <- matrix(T, nrow = length(x_grid), ncol = length(y_grid))

  for (i in 2:(length(x_grid) - 1)) {

    for (j in 2:(length(y_grid) - 1)) {

      isAnyPointInBandRight <- isPointInBandRight(XY_sp, x_grid, y_grid, i - 1, j - 1)

      isAnyPointInBandLeft <- isPointInBandLeft(XY_sp, x_grid, y_grid, i - 1, j - 1)

      isAnyPointInBandUp <- isPointInBandUp(XY_sp, x_grid, y_grid, i - 1, j - 1)

      isAnyPointInBandDown <- isPointInBandDown(XY_sp, x_grid, y_grid, i - 1, j - 1)

      if(!isAnyPointInBandRight | !isAnyPointInBandLeft | !isAnyPointInBandUp | !isAnyPointInBandDown){
        pointInGrid[i,j] <- F
      }

    }

  }

  pointInGrid <- pointInGrid[-c(1,nrow(pointInGrid)),]
  pointInGrid <- pointInGrid[,-c(1,ncol(pointInGrid))]
  x_grid <- x_grid[-c(1,length(x_grid))]
  y_grid <- y_grid[-c(1,length(y_grid))]

  allPoints <- cbind(expand.grid(x_grid, y_grid), as.vector((pointInGrid)))
  allPoints <- allPoints[allPoints[,3],-3]

  allPoints
}

computeSpatialSummaries <- function(Xs, ps, maxPoints){

  n <- nrow(Xs)

  if(ps > 0){

    # isolate unique locations and assign indexes to sites
    uniqueXs <- which(!duplicated(Xs))
    X_s <- Xs[uniqueXs,]

    # indexes assigning original locations (Xs) to new locations (X_s)
    Xs_index <- match(
      do.call(paste, as.data.frame(Xs)),
      do.call(paste, as.data.frame(X_s))
    )

    # indexes assigning the new location X_s to the original
    X_s_index <- match(
      do.call(paste, as.data.frame(X_s)),
      do.call(paste, as.data.frame(Xs))
    )

    # location of support points
    # X_tilde <- as.matrix(buildGrid(X_s, gridStep = .4))

    # assign ps again based on new locations
    if(ps > (nrow(X_s)-1)){
      ps <- nrow(X_s) - 1
    }

    maxPoints <- min(maxPoints, ps)

    list_kmeans <- kmeans(X_s, centers = ps)
    X_tilde <- list_kmeans$centers

    {
      # ggplot() +
      #   geom_point(data = NULL, aes(x = X_tilde[,1],
      #                               y = X_tilde[,2], color = "red")) +
      #   geom_point(data = NULL, aes(x = X_s[,1],
      #                               y = X_s[,2]))

    }

    # distance from support points
    X_s_Xtilde_dist <- t(apply(X_s, 1, function(x){
      apply(X_tilde, 1, function(y){
        (x[1] - y[1])^2 + (x[2] - y[2])^2
      })
    }))

    # indexes of closest support points to the unique points
    X_s_centers <- t(apply(X_s_Xtilde_dist, 1, function(x){
      order(x)[1:maxPoints]
    }))

    # closest support points to the original locations
    Xs_centers <- X_s_centers[Xs_index,]

  } else {

    Xs_index <- rep(NA, nrow(Xs))
    X_s_index <- NULL
    X_tilde <- matrix(NA, ps, 0)
    X_s_centers <- matrix(NA, nrow(Xs), 0)
    Xs_centers <- matrix(NA, nrow(Xs), 0)
    X_s <- matrix(NA, nrow(Xs), 2)

  }

  list(
    "Xs_index" = Xs_index, # indexes matching original locations to new locations
    "X_s_index" = X_s_index, # indexes recovering the new unique locations from the original
    "X_tilde" = X_tilde, # location of support points
    "X_s_centers" = X_s_centers, # indexes to match new locations to the coefficients
    "Xs_centers" = Xs_centers, # indexes to match original locations to the coefficients
    "X_s" = X_s, # new unique locations
    "ps" = ps) # new number of support points

}

precomputeSORmatrices <- function(l_s_grid, list_Xs){

  length_grid_ls <- length(l_s_grid)

  X_tilde <- list_Xs$X_tilde
  X_s <- list_Xs$X_s
  Xs_index <- list_Xs$Xs_index
  X_s_centers <- list_Xs$X_s_centers

  X_centers <- nrow(X_tilde)
  maxPoints <- ncol(X_s_centers)
  n <- length(Xs_index)
  ns <- nrow(X_s)

  Ks_all <- array(NA, dim = c(n, maxPoints, length_grid_ls))
  logDetKuu_grid <- rep(NA, length_grid_ls)
  Lm1_grid <- array(NA, c(ns, ns, length_grid_ls))


  if(X_centers > 0){

    message("Precomputing spatial covariance matrices")

    for (j in 1:length_grid_ls) {

      # if(F){
      #   print(paste0("Precomputing covariance matrix ",j," out of ",length_grid_ls))
      # }

      l_s_current <- l_s_grid[j]

      list_SoRelem <- computeSORmatrix(l_s_current, X_tilde, X_s, Xs_index, X_s_centers)

      Ks_all[,,j] <- list_SoRelem$Ks
      logDetKuu_grid[j] <- list_SoRelem$logDetKuu
      Lm1_grid[,,j] <- list_SoRelem$sq_term

    }
  }

  list("Ks_all" = Ks_all,
       "logDetKuu_grid" = logDetKuu_grid,
       "Lm1_grid" = Lm1_grid,
       "l_s_grid" = l_s_grid)

}

computeSORmatrix <- function(l_s, X_tilde, X_s, Xs_index, X_s_centers){

  ps <- nrow(X_tilde)

  if(ps > 0){

    K_uu <- K2(X_tilde, X_tilde, 1, l_s) + diag(10^(-5), nrow = nrow(X_tilde))
    L_Kmm <- FastGP::rcppeigen_get_chol(K_uu)
    invL_Kmm <- FastGP::rcppeigen_invert_matrix(L_Kmm)
    K_staru <- K2(X_s, X_tilde, 1, l_s)
    KnmLmt <- K_staru %*% t(invL_Kmm)
    Ks <- t(sapply(1:nrow(KnmLmt), function(i){ KnmLmt[i,X_s_centers[i,]]}))
    Ks <- Ks[Xs_index,]

    K_xx <- K2(X_s, X_s, 1, l_s) + diag(exp(-10), nrow = nrow(X_s))
    logDetKuu <- sum(log(FastGP::rcppeigen_get_diag(K_xx))) * 2
    sq_term <- FastGP::rcppeigen_get_chol(FastGP::rcppeigen_invert_matrix(K_xx))

  } else {

    Ks <- matrix(NA, nrow = length(Xs_index), 0)
    logDetKuu <- NULL
    sq_term <- NULL

  }

  list("Ks" = Ks,
       "logDetKuu" = logDetKuu,
       "sq_term" = sq_term)

}

getDefaultSupportPoints <- function(n) min(floor(n * 0.2), n - 1)

# DATA SIMULATION -------

sampleEffects <- function(n){

  sample((-1):1, n, replace = T)

}

createCovMatrix <- function(n, p){

  X <- as.data.frame(matrix(rnorm(n * p, sd = 10), n, p))

  num_cat <- round(0.2 * p)
  cat_cols <- sample(1:p, size = num_cat)

  for (col in cat_cols) {
      X[, col] <- sample(as.character(1:3), size = n, replace = TRUE)
  }

  X

}

simulateData <- function(
    n, S, p, g, gt, d, tau, rnb, ds, ns,
    sigma_b, sigma_bs, sigma_ts, sigma_h, sigma_s, l_s,
    useSpatField, usingSplines, model){

  # simulate data
  {
    speciesNames <- as.character(1:S)

    # Fixed effects matrix
    X <- createCovMatrix(n, p)
    if(p > 0) {

      colnames(X) <- paste("EnvCov", seq_len(p))
      X0 <- X
      list_X <- create_covariates_matrix(
        X,
        spline_vars = usingSplines)
      X <- list_X$X
      p0 <- p
      p <- ncol(X)

    }


    # Traits matrix
    Tr <- matrix(rnorm(S * g), S, g)

    # Spatial locations
    Xs <- matrix(runif(ns * 2), ns, 2)
    if(ns < n){
      Xs <- Xs[sample(1:ns, n, replace = T),]
    }

    ps <- 0
  }

  # params
  {

    # Intercepts
    B0 <- rnorm(S, sd = 1)

    # Traits covariate coefficients
    G <- matrix(sampleEffects(p * g), g, p)

    # Unobserved traits
    A <- matrix(rnorm(gt * S), S, gt)

    # Unobserved response to traits
    C <- matrix(sampleEffects(p * gt), gt, p)
    diag(C) <- 1
    C[lower.tri(C)] <- 0

    # Residual environmental covariates variation
    Bt <- matrix(rnorm(p * S, sd = sigma_b), S, p)

    # Environmental covariates coefficient
    B <- t(computeBtcoef(G, Tr, A, C, Bt))

    # Spatial response to traits
    Gs <- matrix(0, g, ps)

    # Unobserved spatial traits
    As <- matrix(0, S, gt)

    # Unobserved response to spatial traits
    Cs <- matrix(0, gt, ps)
    diag(Cs) <- 1
    Cs[lower.tri(Cs)] <- 0

    # Residual spatial variation
    Bst <- matrix(0, S, ps)

    # Spatial covariates coefficient
    Bs <- t(computeBtcoef(Gs, Tr, As, Cs, Bst))

    #  factor scores
    U <- matrix(rnorm(n * d, sd = sigma_h), n, d)

    # factor loadings
    L <- matrix(sampleEffects(d * S), d, S)
    diag(L) <- 1
    L[lower.tri(L)] <- 0

  }

  # generate spatial field
  {
    if(useSpatField){

      K_mat <- K2(Xs, Xs, sigma_s, l_s) + diag(10^(-5), nrow = ns)
      LU <- chol(K_mat)

      if(ds > 0){

        Lambda <- matrix(rnorm(ds * S), ds, S)
        Sigma <- t(Lambda) %*% Lambda + diag(10^(-5), nrow = S)
        Sigma <- cov2cor(Sigma)
        LV <- chol(Sigma)

        Z <- matrix(rnorm(n * S), n, S)

        spatField <- t(LU) %*% Z %*% LV

      } else { # shared spatial field

        LU <- chol(K_mat)
        Z <- rnorm(n)

        spatField <- matrix(t(LU) %*% Z, n, S, byrow = F)

      }

    } else {
      spatField <- matrix(0, n, S)
    }


    Xs_centers <- matrix(NA, 0, 0)
    Ks <- matrix(NA, 0, 0)

  }

  # simulate observations
  list_eta <- computePsiCoef(
    X, Ks, list_Xs$Xs_centers, Tr,
    B0, G, A, C, Bt,
    Gs, As, Cs, Bst,
    U, L)
  eta <- list_eta$eta
  if(useSpatField) eta <- eta + spatField
  XB <- list_eta$XB
  SE <- list_eta$SE
  UL <- list_eta$UL

  # outcomes
  if(model == "continuous"){

    z <- t(sapply(1:n, function(i){
      sapply(1:S, function(j){
        rnorm(1, eta[i,j], tau[j])
      })
    }))

  } else if(model == "binary"){

    psi <- logistic(eta)

    z <- t(sapply(1:n, function(i){
      sapply(1:S, function(j){
        rbinom(1, 1, psi[i,j])
      })
    }))

  } else if (model == "count"){

    mu <- exp(eta)

    z <- t(sapply(1:n, function(i){
      sapply(1:S, function(j){
        rnbinom(1, size = rnb[j], mu = mu[i,j])
      })
    }))

  }

  varPart <- computeVariancePartitioning_stdevs(XB, SE, UL, model)
  if(useSpatField) {
    varPart <- computeVariancePartitioning_stdevs(XB, SE + spatField, UL, model)
  }

  # save true values
  {
    data <- list(
    "z" = z,
    "X" = X0,
    "Tr" = Tr,
    "Xs" = Xs
    )

    trueParams <- list(
      "B0" = B0,
      "B" = B,
      "Bt" = Bt,
      "G" = G,
      "A" = A,
      "C" = C,
      "Bs" = Bs,
      "Gs" = Gs,
      "As" = As,
      "Cs" = Cs,
      "Bst" = Bst,
      "SE" = SE,
      "spatField" = spatField,
      "U" = U,
      "L" = L,
      "sigma_b" = sigma_b,
      "sigma_bs" = sigma_bs,
      "sigma_s" = sigma_s,
      "l_s" = l_s,
      "eta" = eta,
      "tau" = tau,
      "varPart" = varPart
    )
  }

  list("data" = data,
       "trueParams" = trueParams)
}

# MCMC FUN ---------------------

computeBtcoef <- function(G, Tr, A, C, Btilde){

  Tr %*% G + A %*% C + Btilde

}

computePsiCoef <- function(
    X, Ks, Xs_centers, Tr,
    B0, G, A, C, Bt,
    Gs, As, Cs, Bst,
    U, L){

  ps <- ncol(Bst)

  n <- nrow(X)
  S <- length(B0)

  B0_mat <- matrix(B0, n, S, byrow = T)

  B <- t(computeBtcoef(G, Tr, A, C, Bt))

  XB_prod <- X %*% B

  if(ps > 0){
    Bs <- t(computeBtcoef(Gs, Tr, As, Cs, Bst))
    XsBs_prod <- KsBproduct(Ks, Bs, Xs_centers)
  } else {
    XsBs_prod <- matrix(0, nrow(XB_prod), ncol(XB_prod))
  }

  # PUL <- P %*% H %*% L
  PUL <- U %*% L

  eta <- B0_mat + XB_prod + XsBs_prod + PUL

  list("eta" = eta,
       "XB" = B0_mat + XB_prod,
       "SE" = XsBs_prod,
       "UL" = PUL)

}

# sample residual variance of env. covariates coefficients
sample_sigmab <- function(B, Tr, G, A, C, a_sigmab, b_sigmab){

  p <- nrow(B)
  S <- ncol(B)

  Bt <- t(B) - computeBtcoef(G, Tr, A, C, matrix(0, S, p))

  sumsq <- sum(Bt^2)
  n_samples <- p * S

  sqrt(
    rinvgamma_cpp(a_sigmab + (n_samples / 2), b_sigmab + (sumsq / 2))
  )

}

# sample variance of factor scores (U ~ N(0, sigma_h^2), iid across sites/factors)
sample_sigmah <- function(U, a_sigmah, b_sigmah){

  n <- nrow(U)
  d <- ncol(U)

  sumsq <- sum(U^2)
  n_samples <- n * d

  sqrt(
    rinvgamma_cpp(a_sigmah + (n_samples / 2), b_sigmah + (sumsq / 2))
  )

}

# sample variances of responses
sample_tau <- function(z, eta, a_tau, b_tau){

  n <- nrow(z)
  S <- ncol(z)

  sumsqs <- colSums((z - eta)^2)

  tau <- sapply(1:S, function(s){

    sqrt(
      rinvgamma_cpp(a_tau + (n / 2), b_tau + (sumsqs[s] / 2))
    )
  })

  tau
}

# sample size parameter of responses
sample_rnb <- function(z, eta, tune_sd = 5){

  n <- nrow(z)
  S <- ncol(z)

  mu <- exp(eta)

  rnb <- sapply(1:S, function(s){

    r_current <- rnb[s]

    r_star <- exp(rnorm(1, mean = log(r_current), sd = tune_sd))

    ll_star <- sum(dnbinom(z[,s], size = r_star, mu = mu[,s], log = TRUE))
    ll_current <- sum(dnbinom(z[,s], size = r_current, mu = mu[,s], log = TRUE))

    lp_star <- 0#dgamma(r_star, shape = prior_shape, rate = prior_rate, log = TRUE)
    lp_current <- 0#dgamma(r_current, shape = prior_shape, rate = prior_rate, log = TRUE)

    jacobian <- log(r_star) - log(r_current)

    log_alpha <- (ll_star - ll_current) + (lp_star - lp_current) + jacobian

    # 6. Accept or reject
    if (log(runif(1)) < log_alpha) {
      return(r_star)
    } else {
      return(r_current)
    }

  })

  rnb
}

# sample the fixed effects and the factor loadings
sample_BCsL <- function(k, X, H, G, Tr,
                      A, C, sigma_b,
                      Ks, Xs_centers, Gs, As,
                      Omega) {

  p <- ncol(X)
  d <- ncol(U)
  S <- ncol(Omega)
  ps <- nrow(As)

  M_B <- t(computeBtcoef(G, Tr, A, C, matrix(0, S, p)))

  B <- matrix(0, p, S)
  Cs <- matrix(0, gt, S)
  L <- matrix(0, d, S)

  # no spatial factor loadings if there is no spatial effect
  ps <- ifelse(ps > 0, ps, 0)

  if(p + gts + d > 0){

    for (s in 1:S) {

      k_current <- k[,s]

      XU <- cbind(X, U)

      b_current <- c(M_B[,s], rep(0, ps), rep(0, d))
      B_current <- diag(1, nrow = p + gts + d)
      diag(B_current)[seq_len(p)] <- sigma_b^2

      # if(ps == 0){
      #   BL <- sampleB(XU, B_current, b_current, Omega[,s], k_current)
      # } else {
        invB <- diag(1 / diag(B_current), nrow = p + ps + d)

        BCsL <- sampleB_SoR(XU, B_current, b_current, k_current,
                            Omega[,s], Xs_centers, Ks, ps)
        if(gts > 0){
          Cs[seq_len(gt),s] <- BCsL[p + d + seq_len(gts)]
        }

      # }
      B[seq_len(p),s] <- BCsL[seq_len(p)]
      L[seq_len(d),s] <- BCsL[p + seq_len(d)]


    }

  }

  Bt <- t(B) - Tr %*% G - A %*% C

  Bs <- computeBscoef(Tr, Gs, As, Cs)

  list("B" = B,
       "L" = L,
       "Bt" = Bt,
       "Bs" = Bs,
       "Cs" = Cs)
}

# sample the intercepts, fixed effects, spatial fixed effects and the factor loadings
sample_BBsL <- function(k, X, Tr, U,
                        G, A, C, sigma_b,
                        Gs, As, Cs, sigma_bs,
                        Ks, Xs_centers,
                        Omega, model) {

  p <- ncol(X)
  ps <- ncol(Cs)
  d <- ncol(U)
  S <- ncol(Omega)

  M_B <- t(computeBtcoef(G, Tr, A, C, matrix(0, S, p)))
  M_Bs <- t(computeBtcoef(Gs, Tr, As, Cs, matrix(0, S, ps)))

  B0 <- rep(0, S)
  B <- matrix(0, p, S)
  Bs <- matrix(0, ps, S)
  L <- matrix(0, d, S)

  if(1 + p + ps + d > 0){

    B_current <- diag(1, nrow = 1 + p + d + ps)
    diag(B_current)[1 + seq_len(p)] <- sigma_b^2
    diag(B_current)[1 + p + d + seq_len(ps)] <- sigma_bs^2

    invB_current <- diag(1 / diag(B_current), nrow = 1 + p + d + ps)

    for (s in 1:S) {

      if(model == "continuous"){
        k_current <- k[,s] * Omega[,s]
      } else if(model == "binary"){
        k_current <- k[,s]
      }

      XU <- cbind(1, X, U)

      b_current <- c(0, M_B[,s], rep(0, d), rep(0, ps))

      BBsL <- sampleB_SoR(XU, invB_current, b_current, k_current,
                          Omega[,s], Xs_centers, Ks, ps)

      B0[s] <- BBsL[1]
      B[seq_len(p),s] <- BBsL[1 + seq_len(p)]
      L[seq_len(d),s] <- BBsL[1 + p + seq_len(d)]
      Bs[seq_len(ps),s] <- BBsL[1 + p + d + seq_len(ps)]


    }

  }

  Bt <- t(B) - Tr %*% G - A %*% C

  Bts <- t(Bs) - Tr %*% Gs - As %*% Cs

  list("B" = B,
       "Bt" = Bt,
       "L" = L,
       "Bs" = Bs,
       "Bts" = Bts,
       "B0" = B0)
}

# sample traits (observed and unobserved) response to covariates
sample_GC <- function(B, Tr, A, sigma_b){

  p <- nrow(B)
  S <- ncol(B)
  g <- ncol(Tr)
  gt <- ncol(A)

  G <- matrix(0, g, p)
  C <- matrix(0, gt, p)

  if(p > 0){

    tB <- t(B)

    B_prior <- diag(2, nrow = g + gt)
    b_prior <- rep(0, g + gt)

    for (k in 1:p) {

      B_current <- tB[,k]

      if(g + gt > 0){

        TA <- cbind(Tr, A)

        b_prior <- rep(0, g + gt)
        B_prior <- diag(1, nrow = g + gt)

        GC <- sampleBuniv(TA, B_prior, b_prior, B_current, sigma_b)
        G[seq_len(g),k] <- GC[seq_len(g)]
        C[seq_len(gt),k] <- GC[g + seq_len(gt)]

      }

    }

  }



  Bt <- t(B) - Tr %*% G - A %*% C

  list("G" = G,
       "C" = C,
       "Bt" = Bt)
}

# sample unobserved traits
sample_A <- function(B, C, Tr, G, sigma_b){

  p <- nrow(B)
  gt <- nrow(C)
  S <- ncol(B)

  Btilde <- B - t(Tr %*% G)

  A <- matrix(0, S, gt)

  if(p > 0 & gt > 0){

    B_current <- diag(1, nrow = gt)
    b_current <- rep(0, gt)

    for (s in 1:S) {

      A[s,] <- sampleBuniv(t(C), B_current, b_current, Btilde[,s], sigma_b)

    }

  }

  Bt <- t(B) - Tr %*% G - A %*% C

  list("A" = A,
       "Bt" = Bt)

}

# sample factor scores
sample_U <- function(k, L, X, B, SE, B0, Omega, model){

  d <- nrow(L)
  n <- nrow(k)

  U <- matrix(NA, n, d)

  B0_mat <- matrix(B0, n, S, byrow = T)

  if(d > 0){

    if(model == "continuous"){
      k_new <- k - (B0_mat + X %*% B + SE)
    } else if(model == "binary"){
      k_new <- k - Omega * (B0_mat + X %*% B + SE)
    }

    B_current <- diag(1, nrow = d)
    b_current <- rep(0, d)

    for (i in 1:n) {

      U[i,] <- sampleB(t(L), B_current, b_current, Omega[i,], k_new[i,])

    }
  }

  U

}


loglik_spatialEffect <- function(KsBs_s, Lm1, logdet, sigma_s){

  xsq <- t(KsBs_s) %*% Lm1

  n_p <- nrow(Lm1)

  loglikelihood <- - n_p / 2 * (log(2 * pi) + log(sigma_s^2)) - .5 * logdet -
    (1/2) * (1 / sigma_s^2) * (xsq %*% t(xsq))

  loglikelihood
}

# sample scale parameter of spatial field
sample_ls <- function(idx_ls, SE, list_SoRSummaries,
                      a_l_s, b_l_s, sigma_s){

  if(!is.null(list_SoRSummaries)){

    S <- ncol(SE)

    l_s_grid <- list_SoRSummaries$l_s_grid
    ldet_grid <- list_SoRSummaries$logDetKuu_grid
    Lm1_grid <- list_SoRSummaries$Lm1_grid

    if(idx_ls == 1){
      idx_ls_star <- 2
    } else if(idx_ls == length(l_s_grid)){
      idx_ls_star <- length(l_s_grid) - 1
    } else {
      idx_ls_star <- ifelse(runif(1) < .5, idx_ls - 1, idx_ls + 1)
    }

    # current point
    l_s_current <- l_s_grid[idx_ls]

    loglikelihood_current <- sum(
      sapply(1:S, function(s){
        loglik_spatialEffect(SE[,s], Lm1_grid[,,idx_ls], ldet_grid[idx_ls], sigma_s)
      })
    )

    logPrior_current <- dgamma(l_s_current, a_l_s, b_l_s, log = T)

    logposterior_current <- logPrior_current + loglikelihood_current

    # proposed point
    l_s_star <- l_s_grid[idx_ls_star]

    loglikelihood_star <- sum(
      sapply(1:S, function(s){
        loglik_spatialEffect(SE[,s], Lm1_grid[,,idx_ls_star], ldet_grid[idx_ls_star], sigma_s)
      })
    )

    logPrior_star <- dgamma(l_s_star, a_l_s, b_l_s, log = T)

    logposterior_star <- logPrior_star + loglikelihood_star

    if(runif(1) < exp(logposterior_star - logposterior_current)){
      idx_ls <- idx_ls_star
    }
  }

  idx_ls
}

update_jSDMcoef <- function(list_data,
                            list_params,
                            list_priors,
                            list_Xs,
                            list_SoRSummaries,
                            model,
                            computeVarPart){

  # read data
  {
    z <- list_data$z
    X <- list_data$X
    Tr <- list_data$Tr
  }

  # read params
  {
    B0 <- list_params$B0
    B <- list_params$B
    G <- list_params$G
    A <- list_params$A
    C <- list_params$C
    Bt <- list_params$Bt
    Bs <- list_params$Bs
    Gs <- list_params$Gs
    As <- list_params$As
    Cs <- list_params$Cs
    Bst <- list_params$Bst
    U <- list_params$U
    L <- list_params$L
    sigma_b <- list_params$sigma_b
    sigma_bs <- list_params$sigma_bs
    sigma_h <- list_params$sigma_h
    idx_ls <- list_params$idx_ls
    tau <- list_params$tau

    l_s <- list_SoRSummaries$l_s_grid[idx_ls]
    Ks <- list_SoRSummaries$Ks_all[,,idx_ls]
  }

  # read priors
  {
    a_tau <- list_priors$a_tau
    b_tau <- list_priors$b_tau
    a_sigmab <- list_priors$a_sigmab
    b_sigmab <- list_priors$b_sigmab
    a_sigmabs <- list_priors$a_sigmabs
    b_sigmabs <- list_priors$b_sigmabs
    a_sigmah <- list_priors$a_sigmah
    b_sigmah <- list_priors$b_sigmah
    a_l_s <- list_priors$a_l_s
    b_l_s <- list_priors$b_l_s
  }

  # read state variables
  {
    ps <- nrow(Bs)
  }

  # define transformed target variables
  if(model == "continuous"){
    k <- z
  } else if(model == "binary"){
    k <- z - .5
  } else if(model == "count"){
    # k <- k / w
  }

  # compute linear predictor
  list_psiCoef <- computePsiCoef(
    X, Ks, list_Xs$Xs_centers, Tr,
    B0, G, A, C, Bt,
    Gs, As, Cs, Bst,
    U, L)
  psiCoef <- list_psiCoef$eta
  XB <- list_psiCoef$XB
  SE <- list_psiCoef$SE
  UL <- list_psiCoef$UL

  # sample variance of continuous output
  if(model == "continuous"){
    tau <- sample_tau(z, psiCoef, a_tau, b_tau)
  }

  # sample Omega
  if(model == "continuous"){
    Omega <- matrix(1 / tau^2, nrow(z), ncol(z), byrow = T)
  } else if(model == "binary"){
    Omega <- samplePGvariables(psiCoef)
  }

  # sample fixed effects, spatial trait loadings and factor loadings
  list_BBsL <- sample_BBsL(k, X, Tr, U,
                           G, A, C, sigma_b,
                           Gs, As, Cs, sigma_bs,
                           Ks, list_Xs$Xs_centers,
                           Omega, model)
  B <- list_BBsL$B
  Bt <- list_BBsL$Bt
  Bs <- list_BBsL$Bs
  Bst <- list_BBsL$Bts
  L <- list_BBsL$L
  B0 <- list_BBsL$B0

  # update variance of residuals of environmental covariates
  sigma_b <- sample_sigmab(B, Tr, G, A, C, a_sigmab, b_sigmab)
  if(ps > 0){
    sigma_bs <- sample_sigmab(Bs, Tr, Gs, As, Cs, a_sigmabs, b_sigmabs)
  }

  # sample response to traits (observed and unobsered)
  list_GC <- sample_GC(B, Tr, A, sigma_b)
  G <- list_GC$G
  C <- list_GC$C
  Bt <- list_GC$Bt

  # sample unobserved traits
  list_A <- sample_A(B, C, Tr, G, sigma_b)
  A <- list_A$A
  Bt <- list_A$Bt

  # sample spatial response to traits (observed and unobsered)
  if(ps > 0){
    list_GCs <- sample_GC(Bs, Tr, As, sigma_bs)
    Gs <- list_GCs$G
    Cs <- list_GCs$C
    Bst <- list_GCs$Bt
  }

  # sample unobserved spatial traits
  if(ps > 0){
    list_As <- sample_A(Bs, Cs, Tr, Gs, sigma_bs)
    As <- list_As$A
    Bst <- list_As$Bt
  }

  # sample factor scores
  list_psiCoef <- computePsiCoef(
    X, Ks, list_Xs$Xs_centers, Tr,
    B0, G, A, C, Bt,
    Gs, As, Cs, Bst,
    U, L)
  XB <- list_psiCoef$XB
  SE <- list_psiCoef$SE
  U <- sample_U_cpp(k, L, XB, SE, Omega, sigma_h, model)

  # update variance of factor scores
  sigma_h <- sample_sigmah(U, a_sigmah, b_sigmah)

  # sample spatial field scale
  if(ps > 0){
      idx_ls <- sample_ls(idx_ls, SE,
                          list_SoRSummaries,
                          a_l_s, b_l_s, sigma_s = 1)
      l_s <- list_SoRSummaries$l_s_grid[idx_ls]
      Ks <- list_SoRSummaries$Ks_all[,,idx_ls]
  }

  # output variables
  {
    list_psiCoef <- computePsiCoef(
      X, Ks, list_Xs$Xs_centers, Tr,
      B0, G, A, C, Bt,
      Gs, As, Cs, Bst,
      U, L)
    eta <- list_psiCoef$eta
    XB <- list_psiCoef$XB
    SE <- list_psiCoef$SE
    UL <- list_psiCoef$UL
    if(computeVarPart) {
      variancePartitioning <- computeVariancePartitioning_stdevs(XB, SE, UL, model)
    } else {
      variancePartitioning <- NULL
    }

    if(model == "continuous"){
      eta <- eta
      psi <- NULL
    } else if (model == "binary"){
      psi <- logistic(eta)
    }
  }

  # output params

  list_params <- list(
   "B0" = B0,
   "B" = B,
   "G" = G,
   "A" = A,
   "C" = C,
   "Bt" = Bt,
   "Bs" = Bs,
   "Gs" = Gs,
   "As" = As,
   "Cs" = Cs,
   "Bst" = Bst,
   "U" = U,
   "L" = L,
   "sigma_b" = sigma_b,
   "sigma_bs" = sigma_bs,
   "sigma_h" = sigma_h,
   "idx_ls" = idx_ls,
   "tau" = tau,
   "eta" = eta,
   "psi" = psi,
   "variancePartitioning" = variancePartitioning
  )

  list_params

}

# OUTPUT -------------

plotCoefficient <- function(param_output,
                            covName = NULL,
                            idx_output = NULL){

  if(is.null(covName)){
    stop("No name provided")
  }

  # ncov_psi <- fitModel$infos$ncov_psi
  outputNames <- dimnames(param_output)[[3]]
  covariatesNames <- dimnames(param_output)[[2]]
  S <- dim(param_output)[3]
  idxcov <- which(covariatesNames == covName)

  if(length(idxcov) == 0){
    stop("Covariate name not found. If you are using a categorical covariates,
         the name might have changed to code the level. Use
         colnames(fitModel$Tr) to find the new names")
  }

  if(is.null(idx_output)){
    idx_output <- 1:S
  }

  samples_subset <- matrix(param_output[,idxcov, idx_output],
                           dim(param_output)[1], length(idx_output))

  data_plot <- apply(samples_subset, 2, function(x) {
    quantile(x, probs = c(0.025, 0.975))
  }) %>%
    t %>%
    as.data.frame %>%
    mutate(Output = outputNames[idx_output]) %>%
    mutate(OutputOrder = order(`2.5%`))

  orderOutputs <- order(data_plot$`2.5%`)

  data_plot %>%
    ggplot(aes(x =  factor(Output, level = outputNames[orderOutputs]),
               ymin = `2.5%`,
               ymax = `97.5%`)) + geom_errorbar() +
    ggtitle(covName) +
    theme_bw() +
    theme(
      axis.text = element_text(angle = 90,
                               size = 8),
      plot.title = element_text(hjust = .5,
                                size = 15)
    ) + geom_hline(aes(yintercept = 0), color = "red")

}

computePredictiveProbs <- function(jsdm_output,
                                   X_new,
                                   meanX, sdX,
                                   Xs_new,
                                   meanXs, sdXs,
                                   summarised = F,
                                   confidence = .95,
                                   model){

  X_new <- as.matrix(X_new)

    S <- #fitModel$infos$S
    # speciesNames <- fitModel$infos$speciesNames
    n <- nrow(X_new)

    if(is.null(X_new)) {
      X_psi <- fitModel$X_psi
    }

    B0_output <- jsdm_output$B0_output
    B_output <- jsdm_output$B_output
    Bs_output <- jsdm_output$Bs_output
    L_output <- jsdm_output$L_output

    niter <- dim(B_output)[3]
    nchain <- dim(B_output)[4]

    # transform back coefficients to original scale
    {
      list_BB0_output <- transformCoefficients(B0_output, B_output, meanX, sdX)
      B0_output <- list_BB0_output$B0_output
      B_output <- list_BB0_output$B_output
    }

    if(!summarised){

      eta_output <- array(NA, dim = c(nchain * niter, n, S))
      for (chain in 1:nchain) {
        for (iter in 1:niter) {
          psi_output[iter + (chain - 1)*niter,,] <-
            logistic(
              computePsiCoef(
                X_new, Ks_new, Xs_centers_new, Tr,

                beta_psi_output[,,iter,chain], X_ord,
                          beta_ord_output[,,iter,chain],
                          LL_output[,,iter,chain])
            )
        }
      }

    } else {

      conflevels <- c((1 - confidence)/2, .5, (1 + confidence)/2)

      beta_ord_output <- aperm(apply(beta_ord_output, c(1,2), c), c(2,3,1))
      beta_psi_output <- aperm(apply(beta_psi_output, c(1,2), c), c(2,3,1))
      LL_output <- aperm(apply(LL_output, c(1,2), c), c(2,3,1))

      # niter <- dim(beta_ord_output)[3]

      psi_output <- computePsiOutput(
        X_psi,
        beta_psi_output,
        X_ord,
        beta_ord_output,
        LL_output,
        conflevels)

      # psi_output <- array(NA, dim = c(3, n, S))
      # for (i in 1:n) {
      #   for (j in 1:S) {
      #     mcmc_output <- rep(NA, niter)
      #     for (iter in 1:niter) {
      #       mcmc_output[iter] <- logistic(
      #         computePsiE(X_psi[i,,drop=F], beta_psi_output[,j,iter],
      #                     X_ord[i,,drop=F],
      #                     beta_ord_output[,,iter],
      #                     LL_output[,j,iter])
      #       )
      #
      #     }
      #
      #     psi_output[,i,j] <- quantile(mcmc_output, conflevels)
      #
      #   }
      # }

    }

    psi_output

    # for (iter in 1:niter) {
    #   psi_output[iter,,] <-
    #     logistic(
    #       matrix(beta0_psi_output[iter,], n, S, byrow = T) +
    #         X_psi %*% matrix(beta_psi_output[iter,], ncov_psi, S) +
    #         matrix(U_output[iter,], n, n_factors, byrow = F) %*% matrix(L_output[iter,], n_factors, S)
    #     )
    # }
    #
    # psi_output




}


# valid for logistic regression
pseudo_R2 <- function(eta, y) {

  S <- ncol(y)
  n <- nrow(y)

  mean_y <- mean(y, na.rm = T)

  L_null <- sum(
    sapply(1:S, function(s){
      sapply(1:n, function(i){
        dbinom(y[i,s], 1, mean_y, log = T)
      })
    })
  )

  L_full <- sum(
    sapply(1:S, function(s){
      sapply(1:n, function(i){
        dbinom(y[i,s], 1, logistic(eta[i,s]), log = T)
      })
    })
  )

  1 - L_null / L_full

}

partition_r2 <- function(y, A, B, C) {

  null_ll <- logLik(glm(y ~ 1, family = binomial))

  # Get R² for all subsets
  r2_0 <- 0
  r2_A <- pseudo_R2(A + B + C)
  r2_B <- 1 - logLik(glm(y ~ B, family = binomial)) / null_ll
  r2_C <- 1 - logLik(glm(y ~ C, family = binomial)) / null_ll
  r2_AB <- 1 - logLik(glm(y ~ A + B, family = binomial)) / null_ll
  r2_AC <- 1 - logLik(glm(y ~ A + C, family = binomial)) / null_ll
  r2_BC <- 1 - logLik(glm(y ~ B + C, family = binomial)) / null_ll
  r2_ABC <- 1 - logLik(glm(y ~ A + B + C, family = binomial)) / null_ll

  # LMG averaging (average over all 6 possible orders)
  # A's contribution
  contrib_A <- (r2_A + (r2_AB - r2_B) + (r2_AC - r2_C) + (r2_ABC - r2_BC)) / 4
  contrib_B <- (r2_B + (r2_AB - r2_A) + (r2_BC - r2_C) + (r2_ABC - r2_AC)) / 4
  contrib_C <- (r2_C + (r2_AC - r2_A) + (r2_BC - r2_B) + (r2_ABC - r2_AB)) / 4

  # Convert to percentages (all positive because R² only increases when adding variables)
  total <- contrib_A + contrib_B + contrib_C
  return(c(A = contrib_A/total*100,
           B = contrib_B/total*100,
           C = contrib_C/total*100))
}

computeVariancePartitioning_R2 <- function(XB, SE, UL, y, model){

  S <- ncol(XB)

  if(model == "binary"){
    # Variances for every subset
    R20   <- rep(0, S)
    R2E   <- pseudo_R2(XB, y)
    R2S   <- pseudo_R2(SE, y)
    R2F   <- pseudo_R2(UL, y)
    R2ES  <- pseudo_R2(XB + SE, y)
    R2EF  <- pseudo_R2(XB + UL, y)
    R2SF  <- pseudo_R2(SE + UL, y)
    R2ESF <- pseudo_R2(XB + SE + UL, y)

    # Shapley contributions

    CE <-
      (R2E - R20 +
         (R2ES - R2S) +
         (R2EF - R2F) +
         (R2ESF - R2SF)) / 4

    CS <-
      (R2S - R20 +
         (R2ES - R2E) +
         (R2SF - R2F) +
         (R2ESF - R2EF)) / 4

    CF <-
      (R2F - R20 +
         (R2EF - R2E) +
         (R2SF - R2S) +
         (R2ESF - R2ES)) / 4

    Total <- CE + CS + CF
  } else {
    CE <- NA
    CS <- NA
    CF <- NA
    Total <- NA
  }

  out <- data.frame(
    Environmental = CE / Total,
    Spatial = CS / Total,
    Biotic = CF / Total,
    Total = Total
  )

  out

}

fast_sd <- function(X){
  n <- nrow(X)
  result <- sqrt(colSums((t(t(X) - colMeans(X)))^2) / (n - 1))
  result
}

SD <- function(eta, model) {

  if(model == "continuous"){
    # apply(eta, 2, sd)
    fast_sd(eta)
  } else if (model == "binary"){
    # apply(logistic(eta), 2, sd)
    fast_sd(logistic(eta))
  } else if (model == "count"){
    # apply(exp(eta), 2, sd)
    fast_sd(exp(eta))
  } else {
    stop("Using unavailable model")
  }

}

computeVariancePartitioning_stdevs <- function(XB, SE, UL, model){

  VE   <- SD(XB, model)
  VS   <- SD(SE, model)
  VF   <- SD(UL, model)
  VES  <- SD(XB + SE, model)
  VEF  <- SD(XB + UL, model)
  VSF  <- SD(SE + UL, model)
  VESF <- SD(XB + SE + UL, model)

  VE <- pmax(VE, 0)
  VES_VS <- pmax(VES - VS, 0)
  VEF_VF <- pmax(VEF - VF, 0)
  VESF_VSF <- pmax(VESF - VSF, 0)

  VS <- pmax(VS, 0)
  VES_VE <- pmax(VES - VE, 0)
  VSF_VF <- pmax(VSF - VF, 0)
  VESF_VEF <- pmax(VESF - VEF, 0)

  VF <- pmax(VF, 0)
  VEF_VE <- pmax(VEF - VE, 0)
  VSF_VS <- pmax(VSF - VS, 0)
  VESF_VES <- pmax(VESF - VES, 0)

  CE <- VE + VES_VS + VEF_VF + VESF_VSF

  CS <- VS + VES_VE + VSF_VF + VESF_VEF

  CF <- VF + VEF_VE + VSF_VS + VESF_VES

  Total <- CE + CS + CF

  out <- data.frame(
    Environmental = CE / Total,
    Spatial = CS / Total,
    Biotic = CF / Total,
    Total = VESF
  )

  out

}

computeVariancePartitioning_linvars <- function(XB, SE, UL){

  S <- ncol(XB)

  out <- data.frame(
    Environmental = rep(NA, S),
    Spatial = rep(NA, S),
    Biotic = rep(NA, S),
    Total = rep(NA, S)
  )

  for (s in 1:S) {

    alpha_12 <- var(XB[,s]) / (var(XB[,s]) + var(SE[,s]))
    alpha_23 <- var(SE[,s]) / (var(UL[,s]) + var(SE[,s]))
    alpha_13 <- var(XB[,s]) / (var(XB[,s]) + var(UL[,s]))

    total_var <- var(XB[,s]+SE[,s]+UL[,s])
    contrib <- c(var(XB[,s])+ alpha_12 * cov(XB[,s],SE[,s])+ alpha_13 * cov(XB[,s],UL[,s]),
                 var(SE[,s])+ (1 - alpha_12) * cov(XB[,s],SE[,s])+ alpha_23 * cov(SE[,s],UL[,s]),
                 var(UL[,s])+ (1 - alpha_13) *cov(XB[,s],UL[,s])+ (1 - alpha_23) * cov(SE[,s],UL[,s]))
    percent <- contrib / total_var

    out$Environmental[s] <- percent[1]
    out$Spatial[s] <- percent[2]
    out$Biotic[s] <- percent[3]
    out$Total[s] <- total_var

  }

  out

}

returnVariancePartitioningMatrix <-  function(varPart_output, speciesNames){

  varPart_output_vec <- apply(varPart_output, c(1,2), c)
  varPart_output_mean <- apply(varPart_output_vec, c(2,3), mean)
  varPart_output_sd <- apply(varPart_output_vec, 2, function(x){
    sqrt(sum(apply(x[,1:3], 2, function(y){
      var(y)
    })))
  })

  vp <- data.frame(
    Species = speciesNames,
    Env = varPart_output_mean[,1],
    Spatial = varPart_output_mean[,2],
    Biotic = varPart_output_mean[,3],
    StDev = varPart_output_sd
  )

  vp

}

plotVarPart <- function(varPart_output, speciesNames){

  vp <- returnVariancePartitioningMatrix(varPart_output, speciesNames)

  ggtern(vp,
         aes(x = Env,
             y = Biotic,
             z = Spatial)) +
    geom_point(alpha = 0.7) +
    # labs(
    #   L = "Environment",
    #   T = "Biotic",
    #   R = "Spatial"
    # ) +
    theme_bw() +
    theme_showarrows() +
    theme(
      legend.position = "none",
      tern.axis.title.T = element_text(size = 13),
      tern.axis.title.L = element_text(size = 13),
      tern.axis.title.R = element_text(size = 13)
    )

}

plotCovariateTrend <- function(idxCov, idxSpecies, X0, B, list_ns){

  covNames <- colnames(X0)

  gridVals <- seq(
    min(X0[,idxCov]),
    max(X0[,idxCov]),
    length.out = 100)

  X_new <- createSplinesMatrixSingleCov(list_ns[[idxCov]],
                                        matrix(gridVals, length(gridVals), 1))

  df <- ncol(X_new)

  idxCoeffs <- (idxCov - 1) * df + 1:df
  covEffect <- as.data.frame(as.matrix(X_new) %*% B[idxCoeffs, idxSpecies] )
  covEffect$x <- gridVals

  df_long <- pivot_longer(
    covEffect,
    cols = -x,
    names_to = "Species",
    values_to = "Value"
  )

  ggplot(df_long, aes(x = x, y = Value, colour = Species)) +
    geom_line() +
    theme_bw() + xlab(covNames[idxCov]) + ylab("Covariate Effect")
}

returnSpatialEffectMean <- function(Bs_output, Ks){

  Bs_output_vec <- apply(Bs_output, c(1,2), c)
  Bs_output_vec <- aperm(Bs_output_vec, c(2,3,1))

  spatEffect_mean <- spatEffectMeanCpp(Bs_output_vec, Ks, Xs_centers)

  spatEffect_mean

}

plotSpatialEffect <- function(spatEffect_output, Xs, idx_species = 1){

  spatEffect_output <- returnSpatialEffect(Bs_output, Ks)

  spatEffect_median <- apply(spatEffect_output, c(2,3),
                             function(x) {quantile(x, probs = .5)})

  ggplot(data = NULL, aes(x = Xs[,1],
                          y = Xs[,2],
                          color = spatEffect_median[,idx_species])) + geom_point()

}

returnCorrelationMatrixOutput <- function(L_output, idx_species, outputNames){

  d <- dim(L_output)[1]
  S <- dim(L_output)[2]

  if(is.null(idx_species)){
    idx_species <- 1:S
  }

  L_output_vec <- apply(L_output, c(1,2), c)

  niter <- dim(L_output_vec)[1]

  Lambda_output <- array(NA, dim = c(niter, S, S))

  for (iter in 1:niter) {
    L_output_current <- matrix(L_output_vec[iter,,], S, d, byrow = T)

    Lambda_output[iter,,] <- cov2cor(L_output_current %*% t(L_output_current))
  }

  dimnames(Lambda_output)[[2]] <- outputNames
  dimnames(Lambda_output)[[3]] <- outputNames

  Lambda_output[,idx_species, idx_species]

}

plotCorrelationMatrix <- function(L_output,
                                  idx_species,
                                  speciesNames,
                                  showSignificance = T,
                                  confidence = .95){

  Sigma_output <- returnCorrelationMatrixOutput(L_output, idx_species, speciesNames)

  conflevels <- c((1 - confidence)/2, .5, (1 + confidence)/2)

  Sigma_quantiles <- apply(Sigma_output, c(2,3),
                            function(x){quantile(x, probs = conflevels)})

  sig_matrix <- matrix(FALSE,
                       nrow = dim(Sigma_quantiles)[2],
                       ncol = dim(Sigma_quantiles)[3])
  for(i in 1:nrow(sig_matrix)) {
    for(j in 1:ncol(sig_matrix)) {
      if(i < j) { # Lower triangle only
        lower_bound <- Sigma_quantiles[1,i,j]
        upper_bound <- Sigma_quantiles[3,i,j]
        sig_matrix[i,j] <- (lower_bound < 0 & upper_bound > 0)
      }
    }
  }

  sig_coords <- which(sig_matrix, arr.ind = TRUE)

  p <- ggcorrplot::ggcorrplot(Sigma_quantiles[2,,],
                              method = "square",
                              type = "lower",
                         lab = F, lab_size = 3,
                         colors = c("blue", "white", "red"),
                         title = "Covariance Matrix (as Correlation)") +
    theme(plot.title = element_text(hjust = 0.5,
                                    size = 16,
                                    face = "bold"))

  if(showSignificance){

    annotation_data <- sig_coords[sig_coords[, 1] < sig_coords[, 2], , drop = FALSE]

    # Convert to a data frame for ggplot
    plot_annotations <- data.frame(
      x = annotation_data[, 2] - 1,   # j - 1
      y = annotation_data[, 1],       # i
      label = "x"
    )

    # 2. Add everything to your plot in ONE single layer
    p <- p + geom_text(
      data = plot_annotations,
      aes(x = x, y = y, label = label),
      size = 6,
      color = "black",
      fontface = "bold",
      inherit.aes = FALSE
    )

    # for(k in 1:nrow(sig_coords)) {
    #   i <- sig_coords[k, 1]
    #   j <- sig_coords[k, 2]
    #   # Only add if in lower triangle and i != j
    #   if(i < j) {
    #     p <- p + annotate("text",
    #                       x = j-1,
    #                       y = i,
    #                       label = "x",
    #                       size = 6,
    #                       color = "black",
    #                       fontface = "bold")
    #   }
    # }

  }

  p
}

returnFactorScores <- function(jsdm_output,
                               confidence = .95){

  U_output <- jsdm_output$U_output
  U_output_vec <- apply(U_output, c(1,2), c)

  conflevels <- c((1 - confidence)/2, .5, (1 + confidence)/2)

  U_output_quantiles <- apply(U_output_vec, c(2,3),
                           function(x){quantile(x, probs = conflevels)})

  U_output_quantiles

}

plotFactorScores <- function(jsdm_output,
                             idx_factors,
                             siteNames,
                             confidence = .95){

  factorScoresOutput <- returnFactorScores(jsdm_output, confidence)

  if(length(idx_factors) != 2) {
    stop("Two factors can be used in this plot")
  }

  factorScoresOutput <- factorScoresOutput[,,idx_factors]

  x_coords <- factorScoresOutput[2, , 1]
  y_coords <- factorScoresOutput[2, , 2]
  x_lower <- factorScoresOutput[1, , 1]
  y_lower <- factorScoresOutput[1, , 2]
  x_upper <- factorScoresOutput[3, , 1]
  y_upper <- factorScoresOutput[3, , 2]

  variability <- sqrt((x_upper - x_lower)^2 + (y_upper - y_lower)^2)

  w_x <- x_upper - x_lower
  w_y <- y_upper - y_lower

  plot_data <- data.frame(
    x = x_coords,
    y = y_coords,
    variability = variability,
    r = sqrt((w_x * w_y) / pi),
    name = siteNames
  )

  p <- ggplot(plot_data, aes(x = x, y = y)) +
    # geom_point(aes(size = variability), alpha = 0.6, color = "darkblue") +
    ggforce::geom_circle(aes(x0 = x, y0 = y, r = r),
                fill = "darkblue", alpha = 0.2, color = "darkblue") +
    # geom_text_repel(
    geom_text(
      aes(label = name),
      size = 3.5            # Size of the text font
      # max.overlaps = 20,     # Adjusts how many overlaps it tolerates before hiding some
      # box.padding = 0.3      # Space around the text box
    ) +
    scale_size_area(max_size = 8, name = "Standard Deviation") +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(
      x = paste0("Factor ",idx_factors[1]),
      y = paste0("Factor ",idx_factors[2])
    )

  p
}

returnFactorLoadings_jsdm <- function(jsdm_output,
                               confidence = .95){

  L_output <- jsdm_output$L_output
  L_output_vec <- apply(L_output, c(1,2), c)

  conflevels <- c((1 - confidence)/2, .5, (1 + confidence)/2)

  L_output_quantiles <- apply(L_output_vec, c(2,3),
                           function(x){quantile(x, probs = conflevels)})

  L_output_quantiles

}

plotFactorLoadings_jsdm <- function(jsdm_output,
                                    idx_factors,
                                    speciesNames,
                                    confidence = .95){

  factorLoadingsOutput <- returnFactorLoadings_jsdm(jsdm_output, confidence)

  if(length(idx_factors) != 2) {
    stop("Two factors can be used in this plot")
  }

  factorLoadingsOutput <- factorLoadingsOutput[,idx_factors,]

  x_coords <- factorLoadingsOutput[2, 1, ]
  y_coords <- factorLoadingsOutput[2, 2, ]
  x_lower <- factorLoadingsOutput[1, 1, ]
  y_lower <- factorLoadingsOutput[1, 2, ]
  x_upper <- factorLoadingsOutput[3, 1, ]
  y_upper <- factorLoadingsOutput[3, 2, ]

  variability <- sqrt((x_upper - x_lower)^2 + (y_upper - y_lower)^2)

  w_x <- x_upper - x_lower
  w_y <- y_upper - y_lower

  plot_data <- data.frame(
    x = x_coords,
    y = y_coords,
    variability = variability,
    r = sqrt((w_x * w_y) / pi),
    name = speciesNames
  )

  p <- ggplot(plot_data, aes(x = x, y = y)) +
    # geom_point(aes(size = variability), alpha = 0.6, color = "darkblue") +
    ggforce::geom_circle(aes(x0 = x, y0 = y, r = r),
                fill = "darkblue", alpha = 0.2, color = "darkblue") +
    # geom_text_repel(
    geom_text(
      aes(label = name),
      size = 3.5            # Size of the text font
      # max.overlaps = 20,     # Adjusts how many overlaps it tolerates before hiding some
      # box.padding = 0.3      # Space around the text box
    ) +
    scale_size_area(max_size = 8, name = "Standard Deviation") +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(
      x = paste0("Factor ",idx_factors[1]),
      y = paste0("Factor ",idx_factors[2])
    )

  p
}

# DEPRECATED ---------

# sample traits (observed and unobserved) response to covariates
sample_GC_fixed <- function(B, Tr, A, sigma_b){

  p <- nrow(B)
  S <- ncol(B)
  g <- ncol(Tr)
  gt <- ncol(A)

  tB <- t(B)

  G <- matrix(0, g, p)
  C <- matrix(0, gt, p)

  B_prior <- diag(2, nrow = g + gt)
  b_prior <- rep(0, g + gt)

  for (k in 1:p) {

    elemZero <- max(gt - k, 0)
    elem1 <- ifelse(k <= gt, 1, 0)
    elemNonZero <- gt - elemZero - elem1

    if(elem1 > 0){

      B_current <- tB[,k] - A[,k]

    } else {

      B_current <- tB[,k]

    }

    if(g + elemNonZero > 0){

      TA <- cbind(Tr,
                  matrix(A[,seq_len(elemNonZero)], S, elemNonZero))

      b_prior <- rep(0, g + elemNonZero)
      B_prior <- diag(1, nrow = g + elemNonZero)

      GC <- sampleBuniv(TA, B_prior, b_prior, B_current, sigma_b)
      G[seq_len(g),k] <- GC[seq_len(g)]
      C[seq_len(elemNonZero),k] <- GC[g + seq_len(elemNonZero)]

    }

    if(elem1 > 0) C[k,k] <- 1

  }

  Bt <- t(B) - Tr %*% G - A %*% C

  list("G" = G,
       "C" = C,
       "Bt" = Bt)
}

# multivariate sample of a matrix normal regression
sampleB_m <- function(k, X, eta, Omega, B, b){

  S <- ncol(Omega)
  B_output <- matrix(NA, length(b), S)

  for (s in 1:S) {

    k_current <- k[,s] - Omega[,s] * eta[,s]

    B_output[,s] <- sampleB(X, B, b, Omega[,s], k_current)

  }

  B_output
}

# sample the ordination covariate coefficients
sample_Br <- function(k, L, Xr, Utilde, eta, Omega){

  ncov_ord <- ncol(Xr)
  d <- nrow(L)
  S <- ncol(L)

  Omega_vec <- as.vector(Omega)

  ktilde <- k - (eta + Utilde %*% L) * Omega

  X <- kronecker(t(L), Xr)

  if(ncov_ord > 0 & d > 0){
    Br <- sample_beta_cpp(X,
                          diag(100, nrow = ncol(X)),
                          rep(0, ncol(X)),
                          Omega_vec,
                          ktilde)

  }

  matrix(Br, ncov_ord, d, byrow = F)
}


# sample the fixed effects and the factor loadings
sample_BL_fixed <- function(k, X, H, G, Tr,
                            A, C, sigma_b,
                            Omega) {

  p <- ncol(X)
  d <- ncol(U)
  S <- ncol(Omega)

  M_B <- t(computeBtcoef(G, Tr, A, C, matrix(0, S, p)))

  B <- matrix(0, p, S)
  L <- matrix(0, d, S)

  for (s in 1:S) {

    # elemZero <- max(d - s, 0)
    # elem1 <- ifelse(s <= d, 1, 0)
    # elemNonZero <- d - elemZero - elem1
    #
    # if(elem1 > 0){
    #
    #   k_current <- k[,s] - Omega[,s] * U[,s]
    #
    # } else {
    #
    k_current <- k[,s]
    #
    # }

    elemNonZero <- d

    if(p + elemNonZero > 0){

      XU <- cbind(X, matrix(U[,seq_len(elemNonZero)], n, elemNonZero))

      b_current <- c(M_B[,s], rep(0, elemNonZero))
      B_current <- diag(1, nrow = p + elemNonZero)
      diag(B_current)[seq_len(p)] <- sigma_b^2

      BL <- sampleB(XU, B_current, b_current, Omega[,s], k_current)
      B[seq_len(p),s] <- BL[seq_len(p)]
      L[seq_len(elemNonZero),s] <- BL[p + seq_len(elemNonZero)]

    }

    # if(elem1 > 0) L[s,s] <- 1

  }

  Bt <- t(B) - Tr %*% G - A %*% C

  list("B" = B,
       "L" = L,
       "Bt" = Bt)
}

