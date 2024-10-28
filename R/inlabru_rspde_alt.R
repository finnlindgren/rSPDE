#' @noRd
# Function to process bru's formula

process_formula <- function(bru_result) {
  form <- bru_result$bru_info$model$formula[3]
  form <- as.character(form)
  form <- strsplit(form, "f\\(")
  form <- form[[1]]
  form <- form[-1]
  form_proc <- sub(",.*", "", strsplit(form, "f\\(")[1])
  if (length(form) > 1) {
    for (i in 2:(length(form))) {
      form_proc <- paste(form_proc, " + ", sub(",.*", "", strsplit(form, "f\\(")[i]))
    }
  }
  form_proc <- paste("~", "linkfuninv(", form_proc, ")")
  return(stats::as.formula(form_proc))
}

#' @noRd
# Function to process the link function

process_link <- function(link_name) {
  return_link <- switch(link_name,
    "log" = function(x) {
      INLA::inla.link.log(x, inverse = TRUE)
    },
    "invlog" = function(x) {
      INLA::inla.link.invlog(x, inverse = TRUE)
    },
    "logit" = function(x) {
      INLA::inla.link.logit(x, inverse = TRUE)
    },
    "invlogit" = function(x) {
      INLA::inla.link.invlogit(x, inverse = TRUE)
    },
    "probit" = function(x) {
      INLA::inla.link.probit(x, inverse = TRUE)
    },
    "invprobit" = function(x) {
      INLA::inla.link.invprobit(x, inverse = TRUE)
    },
    "cloglog" = function(x) {
      INLA::inla.link.cloglog(x, inverse = TRUE)
    },
    "invcloglog" = function(x) {
      INLA::inla.link.invcloglog(x, inverse = TRUE)
    },
    "tan" = function(x) {
      INLA::inla.link.tan(x, inverse = TRUE)
    },
    "invtan" = function(x) {
      INLA::inla.link.invtan(x, inverse = TRUE)
    },
    "identity" = function(x) {
      INLA::inla.link.identity(x, inverse = TRUE)
    },
    "invidentity" = function(x) {
      INLA::inla.link.invidentity(x, inverse = TRUE)
    }
  )
  return(return_link)
}

#' @noRd

bru_rerun_with_data <- function(result, idx_data, true_CV, fit_verbose) {
  stopifnot(inherits(result, "bru"))
  if (!true_CV) {
    options <- list(control.mode = list(
      theta = result$mode$theta,
      fixed = TRUE
    ))
  } else {
    options <- list()
  }

  if (fit_verbose) {
    options$verbose <- TRUE
  } else {
    options$verbose <- FALSE
  }

  info <- result[["bru_info"]]
  info[["options"]] <- inlabru::bru_call_options(
    inlabru::bru_options(
      info[["options"]],
      inlabru::as.bru_options(options)
    )
  )

  original_timings <- result[["bru_timings"]]

  lhoods_tmp <- info[["lhoods"]]
  lhoods_tmp[[1]]$response_data$BRU_response[-idx_data] <- NA

  result <- inlabru::iinla(
    model = info[["model"]],
    lhoods = lhoods_tmp,
    initial = result,
    options = info[["options"]]
  )

  new_timings <- result[["bru_iinla"]][["timings"]]$Iteration >
    max(original_timings$Iteration)
  result$bru_timings <-
    rbind(
      original_timings,
      result[["bru_iinla"]][["timings"]][new_timings, , drop = FALSE]
    )

  # Add bru information to the result
  result$bru_info <- info
  class(result) <- c("bru", class(result))
  return(result)
}

#' @noRd

get_post_var <- function(density_df) {
  min_x <- min(density_df[, "x"])
  max_x <- max(density_df[, "x"])
  denstemp <- function(x) {
    dens <- sapply(x, function(z) {
      if (z < min_x) {
        return(0)
      } else if (z > max_x) {
        return(0)
      } else {
        return(approx(x = density_df[, "x"], y = density_df[, "y"], xout = z)$y)
      }
    })
    return(dens)
  }

  post_var <- stats::integrate(
    f = function(z) {
      denstemp(z) * 1 / z
    }, lower = min_x, upper = max_x,
    subdivisions = nrow(density_df),
    stop.on.error = FALSE
  )$value

  return(post_var)
}

#' @noRd

prepare_df_pred <- function(df_pred, result, idx_test) {
  info <- result[["bru_info"]]
  list_of_components <- names(info[["model"]][["effects"]])
  lhoods_tmp <- info[["lhoods"]]

  for (comp in list_of_components) {
    name_input_group <- info[["model"]][["effects"]][[comp]][["group"]][["input"]][["input"]]
    if (!is.null(name_input_group)) {
      name_input_group <- as.character(name_input_group)
      comp_group_tmp <- info[["model"]][["effects"]][[comp]][["env"]][[name_input_group]]
      if (!is.null(comp_group_tmp)) {
        if (!is.null(dim(comp_group_tmp))) {
          comp_group_tmp <- comp_group_tmp[idx_test, , drop = FALSE]
        } else {
          comp_group_tmp <- comp_group_tmp[idx_test]
        }
      } else {
        if (!is.null(dim(lhoods_tmp[[1]]$data[[name_input_group]]))) {
          comp_group_tmp <- lhoods_tmp[[1]]$data[[name_input_group]][idx_test, , drop = FALSE]
        } else {
          comp_group_tmp <- lhoods_tmp[[1]]$data[[name_input_group]][idx_test]
        }
      }
      df_pred[[name_input_group]] <- comp_group_tmp
    }
    name_input_repl <- info[["model"]][["effects"]][[comp]][["replicate"]][["input"]][["input"]]
    if (!is.null(name_input_repl)) {
      name_input_repl <- as.character(name_input_repl)
      comp_repl_tmp <- info[["model"]][["effects"]][[comp]][["env"]][[name_input_repl]]
      if (!is.null(comp_repl_tmp)) {
        if (!is.null(dim(comp_repl_tmp))) {
          comp_repl_tmp <- comp_repl_tmp[idx_test, , drop = FALSE]
        } else {
          comp_repl_tmp <- comp_repl_tmp[idx_test]
        }
      } else {
        if (!is.null(dim(lhoods_tmp[[1]]$data[[name_input_repl]]))) {
          comp_repl_tmp <- lhoods_tmp[[1]]$data[[name_input_repl]][idx_test, , drop = FALSE]
        } else {
          comp_repl_tmp <- lhoods_tmp[[1]]$data[[name_input_repl]][idx_test]
        }
      }
      df_pred[[name_input_repl]] <- comp_repl_tmp
    }
  }
  return(df_pred)
}

#' @noRd

calculate_scores <- function(family, test_data, posterior_samples, hyper_samples, n_samples,
                             parallelize_RP, n_cores_RP, orientation_results) {
  scores <- list()

  if (family == "gaussian") {
    # Calculate MSE
    posterior_mean <- rowMeans(posterior_samples)
    mse <- mean((test_data - posterior_mean)^2)

    # Calculate DSS
    precision_mean <- mean(hyper_samples[, "Precision for the Gaussian observations"])
    Expect_post_var <- 1 / precision_mean
    posterior_variance_of_mean <- rowMeans(posterior_samples^2) - (posterior_mean)^2
    post_var <- Expect_post_var + posterior_variance_of_mean
    dss <- mean((test_data - posterior_mean)^2 / post_var + log(post_var))

    # Calculate CRPS and SCRPS
    phi_sample <- hyper_samples[, "Precision for the Gaussian observations"]
    sd_sample <- 1 / sqrt(phi_sample)

    if (parallelize_RP) {
      cl <- parallel::makeCluster(n_cores_RP)
      doParallel::registerDoParallel(cl)
      Y_sample <- foreach::foreach(i = 1:length(test_data), .combine = 'c') %dopar% {
        posterior_samples[i, ] + sd_sample * rnorm(n_samples)
      }
      parallel::stopCluster(cl)
      Y_sample <- matrix(Y_sample, nrow = length(test_data), byrow = TRUE)
    } else {
      Y_sample <- matrix(0, nrow = length(test_data), ncol = n_samples)
      for (i in 1:length(test_data)) {
        Y_sample[i, ] <- posterior_samples[i, ] + sd_sample * rnorm(n_samples)
      }
    }
    E1_tmp <- rowMeans(abs(Y_sample - test_data))
    E2_tmp <- rowMeans(abs(outer(Y_sample, Y_sample, "-")))
    crps <- mean(-E1_tmp + 0.5 * E2_tmp)
    scrps <- mean(-E1_tmp / E2_tmp - 0.5 * log(E2_tmp))

    scores <- list(mse = mse, dss = dss, crps = crps, scrps = scrps)

  } else if (family == "gamma") {
    # Implement gamma family calculations
    # Similar to the gaussian case, but with appropriate distributions
    # ...

  } else if (family == "poisson") {
    # Implement poisson family calculations
    # ...

  } else {
    stop(paste("Family", family, "is not supported in calculate_scores function."))
  }

  # Adjust scores based on orientation_results
  if (orientation_results == "positive") {
    # For scores where lower is better, negate them
    if ("mse" %in% names(scores)) {
      scores$mse <- -scores$mse
    }
    if ("dss" %in% names(scores)) {
      scores$dss <- -scores$dss
    }
  } else {
    # For scores where higher is better, negate them
    if ("crps" %in% names(scores)) {
      scores$crps <- -scores$crps
    }
    if ("scrps" %in% names(scores)) {
      scores$scrps <- -scores$scrps
    }
  }

  return(scores)
}

#' @name cross_validation
#' @title Perform cross-validation on a list of fitted models.
#' @description Obtain several scores for a list of fitted models according
#' to a folding scheme.
#' @param models A fitted model obtained from calling the `bru()` function or a list of models fitted with the `bru()` function.
#' @param model_names A vector containing the names of the models to appear in the returned `data.frame`. If `NULL`, the names will be of the form `Model 1`, `Model 2`, and so on. By default, it will try to obtain the name from the models list.
#' @param scores A vector containing the scores to be computed. The options are "mse", "crps", "scrps" and "dss". By default, all scores are computed.
#' @param cv_type The type of the folding to be carried out. The options are `k-fold` for `k`-fold cross-validation, in which case the parameter `k` should be provided,
#' `loo`, for leave-one-out and `lpo` for leave-percentage-out, in this case, the parameter `percentage` should be given, and also the `number_folds`
#' with the number of folds to be done. The default is `k-fold`.
#' @param k The number of folds to be used in `k`-fold cross-validation. Will only be used if `cv_type` is `k-fold`.
#' @param percentage The percentage (from 1 to 99) of the data to be used to train the model. Will only be used if `cv_type` is `lpo`.
#' @param number_folds Number of folds to be done if `cv_type` is `lpo`.
#' @param n_samples Number of samples to compute the posterior statistics to be used to compute the scores.
#' @param return_scores_folds If `TRUE`, the scores for each fold will also be returned.
#' @param orientation_results character vector. The options are "negative" and "positive". If "negative", the smaller the scores the better. If "positive", the larger the scores the better.
#' @param include_best Should a row indicating which model was the best for each score be included?
#' @param train_test_indexes A list containing two entries `train`, which is a list whose elements are vectors of indexes of the training data, and `test`, which is a list whose elements are vectors of indexes of the test data.
#' Typically this will be returned list obtained by setting the argument `return_train_test` to `TRUE`.
#' @param return_train_test Logical. Should the training and test indexes be returned? If 'TRUE' the train and test indexes will the 'train_test' element of the returned list.
#' @param return_post_samples If `TRUE` the posterior samples will be included in the returned list.
#' @param parallelize_RP Logical. Should the computation of CRPS and SCRPS (and for some cases, DSS) be parallelized?
#' @param n_cores_RP Number of cores to be used if `parallelize_rp` is `TRUE`.
#' @param true_CV Should a `TRUE` cross-validation be performed? If `TRUE` the models will be fitted on the training dataset. If `FALSE`, the parameters will be kept fixed at the ones obtained in the result object.
#' @param save_settings Logical. If `TRUE`, the settings used in the cross-validation will also be returned.
#' @param print Should partial results be printed throughout the computation?
#' @param fit_verbose Should INLA's run during cross-validation be verbose?
#' @return A data.frame with the fitted models and the corresponding scores.
#' @export

cross_validation <- function(models, model_names = NULL, scores = c("mse", "crps", "scrps", "dss"),
                             cv_type = c("k-fold", "loo", "lpo"),
                             k = 5, percentage = 20, number_folds = 10,
                             n_samples = 1000, return_scores_folds = FALSE,
                             orientation_results = c("negative", "positive"),
                             include_best = TRUE,
                             train_test_indexes = NULL,
                             return_train_test = FALSE,
                             return_post_samples = FALSE,
                             parallelize_RP = FALSE, n_cores_RP = parallel::detectCores() - 1,
                             true_CV = TRUE, save_settings = FALSE,
                             print = TRUE,
                             fit_verbose = FALSE) {
  # Input validation and setup code

  orientation_results <- match.arg(orientation_results)
  scores <- intersect(scores, c("mse", "crps", "scrps", "dss"))
  cv_type <- match.arg(cv_type)

  # Validate models and model_names
  if (inherits(models, "bru")) {
    models <- list(models)
  } else {
    if (!is.list(models)) {
      stop("models must be either a result from a bru call or a list of results from bru() calls!")
    }
    if (!all(sapply(models, inherits, "bru"))) {
      stop("All elements in models must be of class 'bru'!")
    }
  }

  if (is.null(model_names)) {
    model_names <- if (!is.null(names(models))) {
      names(models)
    } else {
      paste("Model", seq_along(models))
    }
  } else {
    if (!is.character(model_names)) {
      stop("model_names must be a character vector!")
    }
    if (length(model_names) != length(models)) {
      stop("model_names must have the same length as models!")
    }
  }

  # Prepare data and indices
  data <- models[[1]]$bru_info$lhoods[[1]]$data
  if (is.vector(data)) data <- as.data.frame(data)

  if (is.null(train_test_indexes)) {
    # Create indices based on cv_type
    if (cv_type == "k-fold") {
      folds <- caret::createFolds(seq_len(nrow(data)), k = k)
      train_list <- lapply(folds, function(x) setdiff(seq_len(nrow(data)), x))
      test_list <- folds
    } else if (cv_type == "loo") {
      train_list <- lapply(seq_len(nrow(data)), function(x) setdiff(seq_len(nrow(data)), x))
      test_list <- as.list(seq_len(nrow(data)))
    } else if (cv_type == "lpo") {
      n_test <- floor((100 - percentage) / 100 * nrow(data))
      train_list <- list()
      test_list <- list()
      for (i in 1:number_folds) {
        test_indices <- sample(seq_len(nrow(data)), n_test)
        train_indices <- setdiff(seq_len(nrow(data)), test_indices)
        train_list[[i]] <- train_indices
        test_list[[i]] <- test_indices
      }
    }
  } else {
    train_list <- train_test_indexes$train
    test_list <- train_test_indexes$test
  }

  n_folds <- length(train_list)
  n_models <- length(models)

  # Initialize score matrices
  score_matrices <- list()
  for (score in scores) {
    score_matrices[[score]] <- matrix(NA_real_, nrow = n_folds, ncol = n_models)
    colnames(score_matrices[[score]]) <- model_names
  }

  # Initialize lists for posterior samples
  if (return_post_samples) {
    post_samples <- setNames(vector("list", n_models), model_names)
    hyper_samples <- setNames(vector("list", n_models), model_names)
    for (m in model_names) {
      post_samples[[m]] <- vector("list", n_folds)
      hyper_samples[[m]] <- vector("list", n_folds)
    }
  }

  # Precompute Formulas
  formula_list <- lapply(models, process_formula)

  # Iterate over folds
  for (fold in seq_len(n_folds)) {
    if (print) cat(sprintf("Processing Fold %d/%d\n", fold, n_folds))

    for (m_idx in seq_len(n_models)) {
      model <- models[[m_idx]]
      model_name <- model_names[m_idx]

      if (print) cat(sprintf("  Evaluating Model: %s\n", model_name))

      test_indices <- test_list[[fold]]
      test_data <- model$bru_info$lhoods[[1]]$response_data$BRU_response[test_indices]

      # Assign link function
      link_name <- model$.args$control.family[[1]]$link
      if (link_name == "default") {
        family <- model$.args$family
        if (family %in% c("gaussian", "gamma", "poisson")) {
          linkfuninv <- function(x) exp(x)
        } else {
          stop(paste("Family", family, "is not supported."))
        }
      } else {
        linkfuninv <- process_link(link_name)
      }

      # Assign link function to formula environment
      formula_tmp <- formula_list[[m_idx]]
      env_tmp <- new.env(parent = environment(formula_tmp))
      assign("linkfuninv", linkfuninv, envir = env_tmp)
      environment(formula_tmp) <- env_tmp

      # Adjust n_samples for specific families
      adjusted_n_samples <- n_samples

      # Perform prediction using group_predict
      post_predict <- group_predict(
        models = model,
        model_names = model_name,
        formula = formula_tmp,
        train_indices = train_list[[fold]],
        test_indices = test_list[[fold]],
        n_samples = adjusted_n_samples,
        pseudo_predict = !true_CV,
        return_samples = return_post_samples,
        return_hyper_samples = return_post_samples,
        n_hyper_samples = 1,
        compute_posterior_means = TRUE,
        print = FALSE,
        fit_verbose = fit_verbose
      )

      # Extract necessary components
      posterior_samples_mat <- post_predict$post_samples
      posterior_mean <- post_predict$post_means
      hyper_samples_1 <- post_predict$hyper_samples
      hyper_marginals <- post_predict$hyper_marginals
      hyper_summary <- post_predict$hyper_summary

      # Store posterior samples if required
      if (return_post_samples) {
        post_samples[[model_name]][[fold]] <- posterior_samples_mat
        hyper_samples[[model_name]][[fold]] <- hyper_samples_1
      }

      # Calculate scores
      scores_list <- calculate_scores(
        family = model$.args$family,
        test_data = test_data,
        posterior_samples = posterior_samples_mat,
        hyper_samples = hyper_samples_1,
        n_samples = n_samples,
        parallelize_RP = parallelize_RP,
        n_cores_RP = n_cores_RP,
        orientation_results = orientation_results
      )

      # Assign scores to matrices
      for (score_name in names(scores_list)) {
        if (score_name %in% scores) {
          score_matrices[[score_name]][fold, m_idx] <- scores_list[[score_name]]
          if (print) cat(sprintf("    %s: %.4f\n", toupper(score_name), scores_list[[score_name]]))
        }
      }
    }
  }

  # Compile Results
  result_df <- data.frame(Model = model_names, stringsAsFactors = FALSE)
  for (score_name in scores) {
    result_df[[score_name]] <- colMeans(score_matrices[[score_name]], na.rm = TRUE)
  }

  # Include Best Models
  if (include_best) {
    best_models <- sapply(result_df[, -1, drop = FALSE], function(col) {
      if (orientation_results == "negative") {
        model_names[which.min(col)]
      } else {
        model_names[which.max(col)]
      }
    })
    best_row <- c("Best", best_models)
    result_df <- rbind(result_df, best_row)
    row.names(result_df)[nrow(result_df)] <- ""
  }

  # Save Settings if Required
  if (save_settings) {
    settings_list <- list(
      n_samples = n_samples,
      cv_type = cv_type,
      true_CV = true_CV,
      orientation_results = orientation_results
    )
    if (cv_type == "k-fold") {
      settings_list$k <- k
    } else if (cv_type == "lpo") {
      settings_list$percentage <- percentage
      settings_list$number_folds <- number_folds
    }
  }

  # Prepare Output
  if (return_scores_folds) {
    out <- list(
      scores_df = result_df,
      scores_folds = score_matrices
    )
    if (save_settings) {
      out$settings <- settings_list
    }
    if (return_train_test) {
      out$train_test <- list(train = train_list, test = test_list)
    }
    if (return_post_samples) {
      out$post_samples <- post_samples
      out$hyper_samples <- hyper_samples
    }
  } else {
    out <- result_df
    if (save_settings) {
      out <- list(scores_df = out, settings = settings_list)
    }
    if (return_train_test) {
      if (is.list(out)) {
        out$train_test <- list(train = train_list, test = test_list)
      } else {
        out <- list(scores_df = out, train_test = list(train = train_list, test = test_list))
      }
    }
    if (return_post_samples) {
      if (is.list(out)) {
        out$post_samples <- post_samples
        out$hyper_samples <- hyper_samples
      } else {
        out <- list(scores_df = out, post_samples = post_samples, hyper_samples = hyper_samples)
      }
    }
  }

  return(out)
}

#' @name group_predict
#' @title Perform prediction on a testing set based on a training set
#' @description Compute prediction of a formula-based expression on a testing set based on a training set.
#' @param models A fitted model obtained from calling the `bru()` function or a list of models fitted with the `bru()` function.
#' @param model_names A vector containing the names of the models to appear in the returned `data.frame`. If `NULL`, the names will be of the form `Model 1`, `Model 2`, and so on. By default, it will try to obtain the name from the models list.
#' @param formula A formula where the right hand side defines an R expression to evaluate for each generated sample. If `NULL``, the latent and hyperparameter states are returned as named list elements. See the manual for the `predict` method in the `inlabru` package.
#' @param train_indices A list containing the indices of the observations for the model to be trained, or a numerical vector containing the indices.
#' @param test_indices A list containing the indices of the test data, where the prediction will be done, or a numerical vector containing the indices.
#' @param n_samples Number of samples to compute the posterior statistics to be used to compute the scores.
#' @param pseudo_predict If `TRUE`, the models will NOT be refitted on the training data, and the parameters obtained on the entire dataset will be used. If `FALSE`, the models will be refitted on the training data.
#' @param return_samples Should the posterior samples be returned?
#' @param return_hyper_samples Should samples for the hyperparameters be obtained?
#' @param n_hyper_samples Number of independent samples of the hyper parameters of size `n_samples`.
#' @param compute_posterior_means Should the posterior means be computed from the posterior samples?
#' @param print Should partial results be printed throughout the computation?
#' @param fit_verbose Should INLA's run during the prediction be verbose?
#' @return A data.frame with the fitted models and the corresponding scores.
#' @export

group_predict <- function(models, model_names = NULL, formula = NULL,
                          train_indices, test_indices, n_samples = 1000,
                          pseudo_predict = TRUE,
                          return_samples = FALSE, return_hyper_samples = FALSE,
                          n_hyper_samples = 1,
                          compute_posterior_means = TRUE,
                          print = TRUE, fit_verbose = FALSE) {
  # Validate inputs
  if (!is.list(models)) {
    stop("models should be a list of 'bru' objects.")
  }
  if (inherits(models, "bru")) {
    models <- list(models)
  } else {
    for (i in 1:length(models)) {
      if (!inherits(models[[i]], "bru")) {
        stop("All models must be 'bru' objects.")
      }
    }
  }

  if (is.null(model_names)) {
    model_names <- names(models)
    if (is.null(model_names)) {
      model_names <- paste0("Model_", seq_along(models))
    }
  }

  if (length(models) != length(model_names)) {
    stop("Length of models and model_names must be the same.")
  }

  # Getting the data
  data <- models[[1]]$bru_info$lhoods[[1]]$data
  if (is.vector(data)) {
    data <- as.data.frame(data)
  }

  post_samples <- list()
  post_means <- list()
  hyper_samples <- list()
  hyper_marginals <- list()
  hyper_summary <- list()

  for (model_number in seq_along(models)) {
    model <- models[[model_number]]
    model_name <- model_names[model_number]

    if (print) {
      cat(sprintf("Processing Model: %s\n", model_name))
    }

    # Prepare training and testing data
    df_train <- prepare_df_pred(data, model, train_indices)
    df_pred <- prepare_df_pred(data, model, test_indices)

    # Refit model if necessary
    if (!pseudo_predict) {
      if (print) cat("Refitting model on training data...\n")
      new_model <- bru_rerun_with_data(model, train_indices, true_CV = TRUE, fit_verbose = fit_verbose)
      if (print) cat("Model refitted.\n")
    } else {
      new_model <- model
    }

    # Generate samples
    if (print) cat("Generating samples...\n")
    generated_samples <- inlabru::generate(new_model, newdata = df_pred, formula = formula, n.samples = n_samples)
    if (print) cat("Samples generated.\n")

    # If only one row, replicate it
    if (nrow(generated_samples) == 1) {
      generated_samples <- matrix(rep(generated_samples, nrow(df_pred)), ncol = ncol(generated_samples), byrow = TRUE)
    }

    post_samples[[model_name]] <- generated_samples

    if (compute_posterior_means) {
      post_means[[model_name]] <- rowMeans(generated_samples)
    }

    hyper_marginals[[model_name]] <- new_model$marginals.hyperpar
    hyper_summary[[model_name]] <- new_model$summary.hyperpar

    if (return_hyper_samples) {
      hyper_samples[[model_name]] <- INLA::inla.hyperpar.sample(n_samples, new_model, improve.marginals = TRUE)
    }
  }

  out <- list()
  if (return_samples) {
    out[["post_samples"]] <- post_samples
  }
  if (return_hyper_samples) {
    out[["hyper_samples"]] <- hyper_samples
  }
  out[["hyper_marginals"]] <- hyper_marginals
  out[["hyper_summary"]] <- hyper_summary
  if (compute_posterior_means) {
    out[["post_means"]] <- post_means
  }
  return(out)
}