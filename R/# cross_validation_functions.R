
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
    },
    stop(paste("Link function", link_name, "is not supported."))
  )
  return(return_link)
}

bru_rerun_with_data <- function(result, idx_data, true_CV, fit_verbose) {
  stopifnot(inherits(result, "bru"))
  if (!true_CV) {
    options <- list(control.mode = list(
      theta = result$mode$theta,
      fixed=TRUE
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

# === Calculate Scores Function ===

calculate_scores <- function(family, test_data, posterior_samples, hyper_samples, n_samples, parallelize_RP, n_cores_RP) {
  scores <- list()
  
  if (family == "gaussian") {
    # Calculate MSE
    posterior_mean <- rowMeans(posterior_samples)
    mse <- mean((test_data - posterior_mean)^2)
    
    # Calculate DSS
    precision_mean <- 1 / mean(hyper_samples[, "Precision for the Gaussian observations"])
    posterior_variance_of_mean <- rowMeans(posterior_samples[, 1:n_samples]^2) - (rowMeans(posterior_samples[, 1:n_samples]))^2
    post_var <- precision_mean + posterior_variance_of_mean
    y_mean <- rowMeans(posterior_samples[, (n_samples + 1):(2 * n_samples)])
    dss <- mean((test_data - y_mean)^2 / post_var + log(post_var))
    
    # Calculate CRPS and SCRPS
    phi_sample_1 <- hyper_samples[, "Precision for the Gaussian observations"][1:n_samples]
    phi_sample_2 <- hyper_samples[, "Precision for the Gaussian observations"][(n_samples + 1):(2 * n_samples)]
    sd_sample_1 <- 1 / sqrt(phi_sample_1)
    sd_sample_2 <- 1 / sqrt(phi_sample_2)
    
    if (parallelize_RP) {
      cl <- makeCluster(n_cores_RP)
      registerDoParallel(cl)
      Y1_sample <- foreach(i = 1:length(test_data), .combine = 'c', .packages = "foreach") %dopar% {
        posterior_samples[i, 1:n_samples] + sd_sample_1[i] * rnorm(n_samples)
      }
      Y2_sample <- foreach(i = 1:length(test_data), .combine = 'c', .packages = "foreach") %dopar% {
        posterior_samples[i, (n_samples + 1):(2 * n_samples)] + sd_sample_2[i] * rnorm(n_samples)
      }
      stopCluster(cl)
      Y1_sample <- split(Y1_sample, rep(1:length(test_data), each = n_samples))
      Y2_sample <- split(Y2_sample, rep(1:length(test_data), each = n_samples))
      
      E1_tmp <- foreach(i = 1:length(test_data), .combine = 'c', .packages = "foreach") %dopar% {
        mean(abs(Y1_sample[[i]] - test_data[i]))
      }
      E2_tmp <- foreach(i = 1:length(test_data), .combine = 'c', .packages = "foreach") %dopar% {
        mean(abs(Y1_sample[[i]] - Y2_sample[[i]]))
      }
    } else {
      Y1_sample <- lapply(1:length(test_data), function(i) {
        posterior_samples[i, 1:n_samples] + sd_sample_1[i] * rnorm(n_samples)
      })
      Y2_sample <- lapply(1:length(test_data), function(i) {
        posterior_samples[i, (n_samples + 1):(2 * n_samples)] + sd_sample_2[i] * rnorm(n_samples)
      })
      
      E1_tmp <- sapply(1:length(test_data), function(i) {
        mean(abs(Y1_sample[[i]] - test_data[i]))
      })
      E2_tmp <- sapply(1:length(test_data), function(i) {
        mean(abs(Y1_sample[[i]] - Y2_sample[[i]]))
      })
    }
    
    crps <- mean(-E1_tmp + 0.5 * E2_tmp)
    scrps <- mean(-E1_tmp / E2_tmp - 0.5 * log(E2_tmp))
    
    scores <- list(mse = mse, dss = dss, crps = crps, scrps = scrps)
    
  } else if (family == "gamma") {
    # Calculate MSE
    posterior_mean <- rowMeans(posterior_samples)
    mse <- mean((test_data - posterior_mean)^2)
    
    # Calculate DSS
    phi_sample_1 <- hyper_samples[, "Precision parameter for the Gamma observations"][1:n_samples]
    phi_sample_2 <- hyper_samples[, "Precision parameter for the Gamma observations"][(n_samples + 1):(2 * n_samples)]
    post_var <- rowMeans(posterior_samples[, 1:n_samples]^2) - (rowMeans(posterior_samples[, 1:n_samples]))^2 +
                rowMeans(posterior_samples[, 1:n_samples])^2 / mean(phi_sample_1)
    y_mean <- rowMeans(posterior_samples[, (n_samples + 1):(2 * n_samples)])
    dss <- mean((test_data - y_mean)^2 / post_var + log(post_var))
    
    # Calculate CRPS and SCRPS
    if (parallelize_RP) {
      cl <- makeCluster(n_cores_RP)
      registerDoParallel(cl)
      Y1_sample <- foreach(i = 1:length(test_data), .combine = 'c', .packages = "foreach") %dopar% {
        scale_temp <- posterior_samples[i, 1:n_samples] / phi_sample_1[i]
        rgamma(n_samples, shape = phi_sample_1[i], scale = scale_temp)
      }
      Y2_sample <- foreach(i = 1:length(test_data), .combine = 'c', .packages = "foreach") %dopar% {
        scale_temp <- posterior_samples[i, (n_samples + 1):(2 * n_samples)] / phi_sample_2[i]
        rgamma(n_samples, shape = phi_sample_2[i], scale = scale_temp)
      }
      stopCluster(cl)
      Y1_sample <- split(Y1_sample, rep(1:length(test_data), each = n_samples))
      Y2_sample <- split(Y2_sample, rep(1:length(test_data), each = n_samples))
      
      E1_tmp <- foreach(i = 1:length(test_data), .combine = 'c', .packages = "foreach") %dopar% {
        mean(abs(Y1_sample[[i]] - test_data[i]))
      }
      E2_tmp <- foreach(i = 1:length(test_data), .combine = 'c', .packages = "foreach") %dopar% {
        mean(abs(Y1_sample[[i]] - Y2_sample[[i]]))
      }
    } else {
      Y1_sample <- lapply(1:length(test_data), function(i) {
        scale_temp <- posterior_samples[i, 1:n_samples] / phi_sample_1[i]
        rgamma(n_samples, shape = phi_sample_1[i], scale = scale_temp)
      })
      Y2_sample <- lapply(1:length(test_data), function(i) {
        scale_temp <- posterior_samples[i, (n_samples + 1):(2 * n_samples)] / phi_sample_2[i]
        rgamma(n_samples, shape = phi_sample_2[i], scale = scale_temp)
      })
      
      E1_tmp <- sapply(1:length(test_data), function(i) {
        mean(abs(Y1_sample[[i]] - test_data[i]))
      })
      E2_tmp <- sapply(1:length(test_data), function(i) {
        mean(abs(Y1_sample[[i]] - Y2_sample[[i]]))
      })
    }
    
    crps <- mean(-E1_tmp + 0.5 * E2_tmp)
    scrps <- mean(-E1_tmp / E2_tmp - 0.5 * log(E2_tmp))
    
    scores <- list(mse = mse, dss = dss, crps = crps, scrps = scrps)
    
  } else if (family == "poisson") {
    # Calculate MSE
    posterior_mean <- rowMeans(posterior_samples)
    mse <- mean((test_data - posterior_mean)^2)
    
    # Calculate DSS
    post_var <- rowMeans(posterior_samples[, 1:n_samples]^2) - (rowMeans(posterior_samples[, 1:n_samples]))^2 +
                rowMeans(posterior_samples[, 1:n_samples])
    y_mean <- rowMeans(posterior_samples[, (n_samples + 1):(2 * n_samples)])
    dss <- mean((test_data - y_mean)^2 / post_var + log(post_var))
    
    # Calculate CRPS and SCRPS
    if (parallelize_RP) {
      cl <- makeCluster(n_cores_RP)
      registerDoParallel(cl)
      Y1_sample <- foreach(i = 1:length(test_data), .combine = 'c', .packages = "foreach") %dopar% {
        rpois(n_samples, posterior_samples[i, 1:n_samples])
      }
      Y2_sample <- foreach(i = 1:length(test_data), .combine = 'c', .packages = "foreach") %dopar% {
        rpois(n_samples, posterior_samples[i, (n_samples + 1):(2 * n_samples)])
      }
      stopCluster(cl)
      Y1_sample <- split(Y1_sample, rep(1:length(test_data), each = n_samples))
      Y2_sample <- split(Y2_sample, rep(1:length(test_data), each = n_samples))
      
      E1_tmp <- foreach(i = 1:length(test_data), .combine = 'c', .packages = "foreach") %dopar% {
        mean(abs(Y1_sample[[i]] - test_data[i]))
      }
      E2_tmp <- foreach(i = 1:length(test_data), .combine = 'c', .packages = "foreach") %dopar% {
        mean(abs(Y1_sample[[i]] - Y2_sample[[i]]))
      }
    } else {
      Y1_sample <- lapply(1:length(test_data), function(i) {
        rpois(n_samples, posterior_samples[i, 1:n_samples])
      })
      Y2_sample <- lapply(1:length(test_data), function(i) {
        rpois(n_samples, posterior_samples[i, (n_samples + 1):(2 * n_samples)])
      })
      
      E1_tmp <- sapply(1:length(test_data), function(i) {
        mean(abs(Y1_sample[[i]] - test_data[i]))
      })
      E2_tmp <- sapply(1:length(test_data), function(i) {
        mean(abs(Y1_sample[[i]] - Y2_sample[[i]]))
      })
    }
    
    crps <- mean(-E1_tmp + 0.5 * E2_tmp)
    scrps <- mean(-E1_tmp / E2_tmp - 0.5 * log(E2_tmp))
    
    scores <- list(mse = mse, dss = dss, crps = crps, scrps = scrps)
    
  } else if (family %in% c("stochvol", "stochvolln", "stochvolnig", "stochvolt")) {
    # Initialize variables based on family
    if (family == "stochvol") {
      # Extract phi parameters
      if ("Offset precision for stochvol" %in% colnames(hyper_samples)) {
        phi_sample_1 <- hyper_samples[, "Offset precision for stochvol"][1:n_samples]
        phi_sample_2 <- hyper_samples[, "Offset precision for stochvol"][(n_samples + 1):(2 * n_samples)]
      } else {
        phi_sample_1 <- NA
        phi_sample_2 <- NA
      }
      
      # Calculate MSE
      posterior_mean <- rowMeans(posterior_samples)
      mse <- mean((test_data - posterior_mean)^2)
      
      # Calculate DSS
      y_mean <- rowMeans(posterior_samples[, (n_samples + 1):(2 * n_samples)])
      post_var <- rowMeans(posterior_samples[, 1:n_samples]^2) - (rowMeans(posterior_samples[, 1:n_samples]))^2 + 
                  1 / phi_sample_1
      dss <- mean((test_data - y_mean)^2 / post_var + log(post_var))
      
      # Calculate CRPS and SCRPS
      if (parallelize_RP) {
        cl <- makeCluster(n_cores_RP)
        registerDoParallel(cl)
        Y1_sample <- foreach(i = 1:length(test_data), .combine = 'c', .packages = "foreach") %dopar% {
          if (is.infinite(phi_sample_1[i])) {
            sqrt(posterior_samples[i, 1:n_samples]) * rnorm(n_samples)
          } else {
            sqrt(posterior_samples[i, 1:n_samples] + 1 / phi_sample_1[i]) * rnorm(n_samples)
          }
        }
        Y2_sample <- foreach(i = 1:length(test_data), .combine = 'c', .packages = "foreach") %dopar% {
          if (is.infinite(phi_sample_2[i])) {
            sqrt(posterior_samples[i, (n_samples + 1):(2 * n_samples)]) * rnorm(n_samples)
          } else {
            sqrt(posterior_samples[i, (n_samples + 1):(2 * n_samples)] + 1 / phi_sample_2[i]) * rnorm(n_samples)
          }
        }
        stopCluster(cl)
        Y1_sample <- split(Y1_sample, rep(1:length(test_data), each = n_samples))
        Y2_sample <- split(Y2_sample, rep(1:length(test_data), each = n_samples))
        
        E1_tmp <- foreach(i = 1:length(test_data), .combine = 'c', .packages = "foreach") %dopar% {
          mean(abs(Y1_sample[[i]] - test_data[i]))
        }
        E2_tmp <- foreach(i = 1:length(test_data), .combine = 'c', .packages = "foreach") %dopar% {
          mean(abs(Y1_sample[[i]] - Y2_sample[[i]]))
        }
      } else {
        Y1_sample <- lapply(1:length(test_data), function(i) {
          if (is.infinite(phi_sample_1[i])) {
            sqrt(posterior_samples[i, 1:n_samples]) * rnorm(n_samples)
          } else {
            sqrt(posterior_samples[i, 1:n_samples] + 1 / phi_sample_1[i]) * rnorm(n_samples)
          }
        })
        Y2_sample <- lapply(1:length(test_data), function(i) {
          if (is.infinite(phi_sample_2[i])) {
            sqrt(posterior_samples[i, (n_samples + 1):(2 * n_samples)]) * rnorm(n_samples)
          } else {
            sqrt(posterior_samples[i, (n_samples + 1):(2 * n_samples)] + 1 / phi_sample_2[i]) * rnorm(n_samples)
          }
        })
        
        E1_tmp <- sapply(1:length(test_data), function(i) {
          mean(abs(Y1_sample[[i]] - test_data[i]))
        })
        E2_tmp <- sapply(1:length(test_data), function(i) {
          mean(abs(Y1_sample[[i]] - Y2_sample[[i]]))
        })
      }
      
      crps <- mean(-E1_tmp + 0.5 * E2_tmp)
      scrps <- mean(-E1_tmp / E2_tmp - 0.5 * log(E2_tmp))
      
      scores <- list(mse = mse, dss = dss, crps = crps, scrps = scrps)
      
    } else if (family == "stochvolln") {
      # Extract relevant parameters
      if ("Offset precision for stochvolln" %in% colnames(hyper_samples)) {
        phi_sample_1 <- hyper_samples[, "Offset precision for stochvolln"][1:n_samples]
        phi_sample_2 <- hyper_samples[, "Offset precision for stochvolln"][(n_samples + 1):(2 * n_samples)]
      } else {
        phi_sample_1 <- NA
        phi_sample_2 <- NA
      }
      
      if ("Mean offset for stochvolln" %in% colnames(hyper_samples)) {
        mu_sample_1 <- hyper_samples[, "Mean offset for stochvolln"][1:n_samples]
        mu_sample_2 <- hyper_samples[, "Mean offset for stochvolln"][(n_samples + 1):(2 * n_samples)]
      } else {
        mu_sample_1 <- NA
        mu_sample_2 <- NA
      }
      
      # Calculate MSE
      posterior_mean <- rowMeans(posterior_samples)
      mse <- mean((test_data - posterior_mean)^2)
      
      # Calculate DSS
      y_mean <- rowMeans(posterior_samples[, (n_samples + 1):(2 * n_samples)])
      post_var <- rowMeans(posterior_samples[, 1:n_samples]^2) - (rowMeans(posterior_samples[, 1:n_samples]))^2 +
                  rowMeans(posterior_samples[, 1:n_samples])^2 / mean(phi_sample_1)
      dss <- mean((test_data - y_mean)^2 / post_var + log(post_var))
      
      # Calculate CRPS and SCRPS
      if (parallelize_RP) {
        cl <- makeCluster(n_cores_RP)
        registerDoParallel(cl)
        Y1_sample <- foreach(i = 1:length(test_data), .combine = 'c', .packages = c("foreach", "ngme2")) %dopar% {
          if (is.infinite(phi_sample_1[i])) {
            mu_sample_1[i] + sqrt(1 / phi_sample_1[i]) * rnorm(n_samples)
          } else {
            mu_sample_1[i] + sqrt(1 / phi_sample_1[i]) * rnorm(n_samples)
          }
        }
        Y2_sample <- foreach(i = 1:length(test_data), .combine = 'c', .packages = c("foreach", "ngme2")) %dopar% {
          if (is.infinite(phi_sample_2[i])) {
            mu_sample_2[i] + sqrt(1 / phi_sample_2[i]) * rnorm(n_samples)
          } else {
            mu_sample_2[i] + sqrt(1 / phi_sample_2[i]) * rnorm(n_samples)
          }
        }
        stopCluster(cl)
        Y1_sample <- split(Y1_sample, rep(1:length(test_data), each = n_samples))
        Y2_sample <- split(Y2_sample, rep(1:length(test_data), each = n_samples))
        
        E1_tmp <- foreach(i = 1:length(test_data), .combine = 'c', .packages = "foreach") %dopar% {
          mean(abs(Y1_sample[[i]] - test_data[i]))
        }
        E2_tmp <- foreach(i = 1:length(test_data), .combine = 'c', .packages = "foreach") %dopar% {
          mean(abs(Y1_sample[[i]] - Y2_sample[[i]]))
        }
      } else {
        Y1_sample <- lapply(1:length(test_data), function(i) {
          if (is.infinite(phi_sample_1[i])) {
            mu_sample_1[i] + sqrt(1 / phi_sample_1[i]) * rnorm(n_samples)
          } else {
            mu_sample_1[i] + sqrt(1 / phi_sample_1[i]) * rnorm(n_samples)
          }
        })
        Y2_sample <- lapply(1:length(test_data), function(i) {
          if (is.infinite(phi_sample_2[i])) {
            mu_sample_2[i] + sqrt(1 / phi_sample_2[i]) * rnorm(n_samples)
          } else {
            mu_sample_2[i] + sqrt(1 / phi_sample_2[i]) * rnorm(n_samples)
          }
        })
        
        E1_tmp <- sapply(1:length(test_data), function(i) {
          mean(abs(Y1_sample[[i]] - test_data[i]))
        })
        E2_tmp <- sapply(1:length(test_data), function(i) {
          mean(abs(Y1_sample[[i]] - Y2_sample[[i]]))
        })
      }
      
      crps <- mean(-E1_tmp + 0.5 * E2_tmp)
      scrps <- mean(-E1_tmp / E2_tmp - 0.5 * log(E2_tmp))
      
      scores <- list(mse = mse, dss = dss, crps = crps, scrps = scrps)
      
    } else if (family == "stochvolnig") {
      # Implement handling for stochvolnig
      # Extract relevant parameters
      if ("shape parameter for stochvol-nig" %in% colnames(hyper_samples)) {
        shape_1 <- hyper_samples[, "shape parameter for stochvol-nig"][1:n_samples]
        shape_2 <- hyper_samples[, "shape parameter for stochvol-nig"][(n_samples + 1):(2 * n_samples)]
      } else {
        shape_1 <- NA
        shape_2 <- NA
      }
      
      if ("skewness parameter for stochvol-nig" %in% colnames(hyper_samples)) {
        skewness_1 <- hyper_samples[, "skewness parameter for stochvol-nig"][1:n_samples]
        skewness_2 <- hyper_samples[, "skewness parameter for stochvol-nig"][(n_samples + 1):(2 * n_samples)]
      } else {
        skewness_1 <- NA
        skewness_2 <- NA
      }
            
      # Calculate MSE
      posterior_mean <- rowMeans(posterior_samples)
      mse <- mean((test_data - posterior_mean)^2)
      
      # Calculate DSS
      y_mean <- rowMeans(posterior_samples[, (n_samples + 1):(2 * n_samples)])
      post_var <- rowMeans(posterior_samples[, 1:n_samples]^2) - (rowMeans(posterior_samples[, 1:n_samples]))^2 +
                  rowMeans(posterior_samples[, 1:n_samples])^2 / mean(phi_sample_1)
      dss <- mean((test_data - y_mean)^2 / post_var + log(post_var))
      
      # Calculate CRPS and SCRPS
      if (parallelize_RP) {
        cl <- makeCluster(n_cores_RP)
        registerDoParallel(cl)
        Y1_sample <- foreach(i = 1:length(test_data), .combine = 'c', .packages = c("foreach", "ngme2")) %dopar% {
          sqrt(posterior_samples[i, 1:n_samples] + 1 / phi_sample_1[i]) * ngme2::rnig(n_samples, 
                                                                                       delta = -skewness_1[i]/sqrt(1 + (skewness_1[i]^2 / shape_1[i]^2)),
                                                                                       mu = skewness_1[i],
                                                                                       nu = shape_1[i]^2,
                                                                                       sigma = 1/sqrt(1 + (skewness_1[i]^2 / shape_1[i]^2)))
        }
        Y2_sample <- foreach(i = 1:length(test_data), .combine = 'c', .packages = c("foreach", "ngme2")) %dopar% {
          sqrt(posterior_samples[i, (n_samples + 1):(2 * n_samples)] + 1 / phi_sample_2[i]) * ngme2::rnig(n_samples, 
                                                                                                           delta = -skewness_2[i]/sqrt(1 + (skewness_2[i]^2 / shape_2[i]^2)),
                                                                                                           mu = skewness_2[i],
                                                                                                           nu = shape_2[i]^2,
                                                                                                           sigma = 1/sqrt(1 + (skewness_2[i]^2 / shape_2[i]^2)))
        }
        stopCluster(cl)
        Y1_sample <- split(Y1_sample, rep(1:length(test_data), each = n_samples))
        Y2_sample <- split(Y2_sample, rep(1:length(test_data), each = n_samples))
        
        E1_tmp <- foreach(i = 1:length(test_data), .combine = 'c', .packages = "foreach") %dopar% {
          mean(abs(Y1_sample[[i]] - test_data[i]))
        }
        E2_tmp <- foreach(i = 1:length(test_data), .combine = 'c', .packages = "foreach") %dopar% {
          mean(abs(Y1_sample[[i]] - Y2_sample[[i]]))
        }
      } else {
        Y1_sample <- lapply(1:length(test_data), function(i) {
          sqrt(posterior_samples[i, 1:n_samples] + 1 / phi_sample_1[i]) * ngme2::rnig(n_samples, 
                                                                                       delta = -skewness_1[i]/sqrt(1 + (skewness_1[i]^2 / shape_1[i]^2)),
                                                                                       mu = skewness_1[i],
                                                                                       nu = shape_1[i]^2,
                                                                                       sigma = 1/sqrt(1 + (skewness_1[i]^2 / shape_1[i]^2)))
        })
        Y2_sample <- lapply(1:length(test_data), function(i) {
          sqrt(posterior_samples[i, (n_samples + 1):(2 * n_samples)] + 1 / phi_sample_2[i]) * ngme2::rnig(n_samples, 
                                                                                                           delta = -skewness_2[i]/sqrt(1 + (skewness_2[i]^2 / shape_2[i]^2)),
                                                                                                           mu = skewness_2[i],
                                                                                                           nu = shape_2[i]^2,
                                                                                                           sigma = 1/sqrt(1 + (skewness_2[i]^2 / shape_2[i]^2)))
        })
        
        E1_tmp <- sapply(1:length(test_data), function(i) {
          mean(abs(Y1_sample[[i]] - test_data[i]))
        })
        E2_tmp <- sapply(1:length(test_data), function(i) {
          mean(abs(Y1_sample[[i]] - Y2_sample[[i]]))
        })
      }
      
      crps <- mean(-E1_tmp + 0.5 * E2_tmp)
      scrps <- mean(-E1_tmp / E2_tmp - 0.5 * log(E2_tmp))
      
      scores <- list(mse = mse, dss = dss, crps = crps, scrps = scrps)
      
    } else if (family == "stochvolt") {
      # Extract relevant parameters
      if ("degrees of freedom for stochvol student-t" %in% colnames(hyper_samples)) {
        degree_1 <- hyper_samples[, "degrees of freedom for stochvol student-t"][1:n_samples]
        degree_2 <- hyper_samples[, "degrees of freedom for stochvol student-t"][(n_samples + 1):(2 * n_samples)]
      } else {
        degree_1 <- NA
        degree_2 <- NA
      }
      
      # Calculate MSE
      posterior_mean <- rowMeans(posterior_samples)
      mse <- mean((test_data - posterior_mean)^2)
      
      # Calculate DSS
      y_mean <- rowMeans(posterior_samples[, (n_samples + 1):(2 * n_samples)])
      post_var <- rowMeans(posterior_samples[, 1:n_samples]^2) - (rowMeans(posterior_samples[, 1:n_samples]))^2 +
                  rowMeans(posterior_samples[, 1:n_samples])
      dss <- mean((test_data - y_mean)^2 / post_var + log(post_var))
      
      # Calculate CRPS and SCRPS
      if (parallelize_RP) {
        cl <- makeCluster(n_cores_RP)
        registerDoParallel(cl)
        Y1_sample <- foreach(i = 1:length(test_data), .combine = 'c', .packages = "foreach") %dopar% {
          sqrt(posterior_samples[i, 1:n_samples]) * rt(n_samples, degree_1[i])
        }
        Y2_sample <- foreach(i = 1:length(test_data), .combine = 'c', .packages = "foreach") %dopar% {
          sqrt(posterior_samples[i, (n_samples + 1):(2 * n_samples)]) * rt(n_samples, degree_2[i])
        }
        stopCluster(cl)
        Y1_sample <- split(Y1_sample, rep(1:length(test_data), each = n_samples))
        Y2_sample <- split(Y2_sample, rep(1:length(test_data), each = n_samples))
        
        E1_tmp <- foreach(i = 1:length(test_data), .combine = 'c', .packages = "foreach") %dopar% {
          mean(abs(Y1_sample[[i]] - test_data[i]))
        }
        E2_tmp <- foreach(i = 1:length(test_data), .combine = 'c', .packages = "foreach") %dopar% {
          mean(abs(Y1_sample[[i]] - Y2_sample[[i]]))
        }
      } else {
        Y1_sample <- lapply(1:length(test_data), function(i) {
          sqrt(posterior_samples[i, 1:n_samples]) * rt(n_samples, degree_1[i])
        })
        Y2_sample <- lapply(1:length(test_data), function(i) {
          sqrt(posterior_samples[i, (n_samples + 1):(2 * n_samples)]) * rt(n_samples, degree_2[i])
        })
        
        E1_tmp <- sapply(1:length(test_data), function(i) {
          mean(abs(Y1_sample[[i]] - test_data[i]))
        })
        E2_tmp <- sapply(1:length(test_data), function(i) {
          mean(abs(Y1_sample[[i]] - Y2_sample[[i]]))
        })
      }
      
      crps <- mean(-E1_tmp + 0.5 * E2_tmp)
      scrps <- mean(-E1_tmp / E2_tmp - 0.5 * log(E2_tmp))
      
      scores <- list(mse = mse, dss = dss, crps = crps, scrps = scrps)
      
    } 
    } else {
      stop(paste("Family", family, "is not supported in calculate_scores function."))
    }
    
    return(scores)
}

# === group_predict Function ===

group_predict <- function(model, model_name, formula = NULL,
                          train_indices, test_indices, n_samples = 1000,
                          pseudo_predict = TRUE,
                          return_samples = FALSE, return_hyper_samples = FALSE,
                          n_hyper_samples = 1,
                          compute_posterior_means = TRUE,
                          print = FALSE, fit_verbose = FALSE) {
  # Validate inputs
  if (missing(model)) stop("model must be provided!")
  if (!inherits(model, "bru")) stop("model must be of class 'bru'!")
  
  # Extract data
  data <- model$bru_info$lhoods[[1]]$data
  if (is.vector(data)) data <- as.data.frame(data)
  
  # Select training and testing data
  train_data <- prepare_df_pred(data, model, train_indices)
  test_data_df <- prepare_df_pred(data, model, test_indices)
  
  # Refit model on training data if pseudo_predict is FALSE
  if (!pseudo_predict) {
    if (print) cat("Refitting model on training data...\n")
    new_model <- bru_rerun_with_data(model, train_indices, true_CV = TRUE, fit_verbose = fit_verbose)
    if (print) cat("Model refitted.\n")
  } else {
    new_model <- model
  }
  
  # Generate posterior samples
  if (print) cat("Generating posterior samples...\n")
  post_predict <- inlabru::generate(new_model, newdata = test_data_df, formula = formula, n.samples = n_samples)
  if (print) cat("Posterior samples generated.\n")
  
  # Extract posterior samples and summary
  posterior_samples <- post_predict$post_samples
  posterior_means <- if (compute_posterior_means) rowMeans(posterior_samples) else NULL
  hyper_samples <- if (return_hyper_samples) {
    INLA::inla.hyperpar.sample(n_hyper_samples, new_model, n.samples = n_samples, improve.marginals = TRUE)
  } else {
    NULL
  }
  
  hyper_marginals <- new_model$marginals.hyperpar
  hyper_summary <- new_model$summary.hyperpar
  
  return(list(
    post_samples = posterior_samples,
    post_means = posterior_means,
    hyper_samples = hyper_samples,
    hyper_marginals = hyper_marginals,
    hyper_summary = hyper_summary
  ))
}

# === cross_validation Function ===

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
  # === Input Validation ===
  
  orientation_results <- match.arg(orientation_results)
  scores <- intersect(scores, c("mse", "crps", "scrps", "dss"))
  cv_type <- match.arg(cv_type)
  
  # Validate percentage and number_folds
  if (!is.numeric(percentage) || length(percentage) != 1) {
    stop("percentage must be a single numeric value.")
  }
  if (percentage %% 1 != 0) {
    warning("Non-integer percentage given, it will be rounded to an integer number.")
    percentage <- round(percentage)
  }
  if (percentage <= 0 || percentage >= 100) {
    stop("percentage must be a number between 1 and 99!")
  }
  
  if (!is.numeric(number_folds) || length(number_folds) != 1) {
    stop("number_folds must be a single numeric value.")
  }
  if (number_folds %% 1 != 0) {
    warning("Non-integer number_folds given, it will be rounded to an integer number.")
    number_folds <- round(number_folds)
  }
  if (number_folds <= 0) {
    stop("number_folds must be positive!")
  }
  
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
  
  # === Parallelization Setup ===
  if (parallelize_RP) {
    cl <- makeCluster(n_cores_RP)
    registerDoParallel(cl)
    on.exit({
      stopCluster(cl)
    }, add = TRUE)
  }
  
  # === Extract Data ===
  data <- models[[1]]$bru_info$lhoods[[1]]$data
  if (is.vector(data)) data <- as.data.frame(data)
  
  # === Create Train-Test Indices ===
  if (is.null(train_test_indexes)) {
    train_test_indices <- create_train_test_indices(data, cv_type, k, percentage, number_folds)
    train_list <- train_test_indices$train
    test_list <- train_test_indices$test
  } else {
    if (!is.list(train_test_indexes) || 
        is.null(train_test_indexes$train) || 
        is.null(train_test_indexes$test)) {
      stop("train_test_indexes must be a list containing 'train' and 'test' elements.")
    }
    if (!is.list(train_test_indexes$train) || !is.list(train_test_indexes$test)) {
      stop("'train' and 'test' elements in train_test_indexes must be lists.")
    }
    train_list <- train_test_indexes$train
    test_list <- train_test_indexes$test
  }
  
  n_folds <- length(train_list)
  n_models <- length(models)
  
  # === Initialize Score Matrices ===
  score_matrices <- list()
  if ("dss" %in% scores) {
    score_matrices$dss <- matrix(NA_real_, nrow = n_folds, ncol = n_models)
    colnames(score_matrices$dss) <- model_names
  }
  if ("mse" %in% scores) {
    score_matrices$mse <- matrix(NA_real_, nrow = n_folds, ncol = n_models)
    colnames(score_matrices$mse) <- model_names
  }
  if ("crps" %in% scores) {
    score_matrices$crps <- matrix(NA_real_, nrow = n_folds, ncol = n_models)
    colnames(score_matrices$crps) <- model_names
  }
  if ("scrps" %in% scores) {
    score_matrices$scrps <- matrix(NA_real_, nrow = n_folds, ncol = n_models)
    colnames(score_matrices$scrps) <- model_names
  }
  
  # === Initialize Lists for Posterior Samples ===
  if (return_post_samples) {
    post_samples <- setNames(vector("list", n_models), model_names)
    hyper_samples <- setNames(vector("list", n_models), model_names)
    for (m in model_names) {
      post_samples[[m]] <- vector("list", n_folds)
      hyper_samples[[m]] <- vector("list", n_folds)
    }
  }
  
  # === Precompute Formulas ===
  formula_list <- lapply(models, process_formula)
  
  # === Iterate Over Folds ===
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
        linkfuninv <- switch(family,
                             "gaussian" = identity,
                             "gamma" = exp,
                             "poisson" = exp,
                             "stochvol" = exp,
                             "stochvolln" = exp,
                             "stochvolnig" = exp,
                             "stochvolt" = exp,
                             stop(paste("The family", family, "is not supported yet.")))
      } else {
        linkfuninv <- process_link(link_name)
      }
      
      # Assign link function to formula environment
      formula_tmp <- formula_list[[m_idx]]
      env_tmp <- new.env(parent = environment(formula_tmp))
      assign("linkfuninv", linkfuninv, envir = env_tmp)
      environment(formula_tmp) <- env_tmp
      
      # Adjust n_samples for specific families
      if (model$.args$family %in% c("stochvol", "stochvolln", "stochvolnig", "stochvolt")) {
        adjusted_n_samples <- 2 * n_samples
      } else {
        adjusted_n_samples <- n_samples
      }
      
      # Perform prediction using group_predict
      post_predict <- group_predict(
        model = model,
        model_name = model_name,
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
      
      # === Calculate Scores based on Family ===
      family <- model$.args$family
      
      scores <- calculate_scores(family, test_data, posterior_samples_mat, hyper_samples_1, 
                                 n_samples, parallelize_RP, n_cores_RP)
      
      # Assign scores to matrices
      if ("mse" %in% scores) {
        mse_val <- scores$mse
        if (orientation_results == "positive") mse_val <- -mse_val
        score_matrices$mse[fold, m_idx] <- mse_val
        if (print) cat(sprintf("    MSE: %.4f\n", mse_val))
      }
      if ("dss" %in% scores) {
        dss_val <- scores$dss
        if (orientation_results == "positive") dss_val <- -dss_val
        score_matrices$dss[fold, m_idx] <- dss_val
        if (print) cat(sprintf("    DSS: %.4f\n", dss_val))
      }
      if ("crps" %in% scores) {
        crps_val <- scores$crps
        if (orientation_results == "negative") crps_val <- -crps_val
        score_matrices$crps[fold, m_idx] <- crps_val
        if (print) cat(sprintf("    CRPS: %.4f\n", crps_val))
      }
      if ("scrps" %in% scores) {
        scrps_val <- scores$scrps
        if (orientation_results == "negative") scrps_val <- -scrps_val
        score_matrices$scrps[fold, m_idx] <- scrps_val
        if (print) cat(sprintf("    SCRPS: %.4f\n", scrps_val))
      }
    }
  }
  
  # === Compile Results ===
  result_df <- data.frame(Model = model_names, stringsAsFactors = FALSE)
  
  if ("dss" %in% scores) {
    dss_mean <- colMeans(score_matrices$dss, na.rm = TRUE)
    result_df$dss <- dss_mean
  }
  if ("mse" %in% scores) {
    mse_mean <- colMeans(score_matrices$mse, na.rm = TRUE)
    result_df$mse <- mse_mean
  }
  if ("crps" %in% scores) {
    crps_mean <- colMeans(score_matrices$crps, na.rm = TRUE)
    result_df$crps <- crps_mean
  }
  if ("scrps" %in% scores) {
    scrps_mean <- colMeans(score_matrices$scrps, na.rm = TRUE)
    result_df$scrps <- scrps_mean
  }
  
  # === Include Best Models ===
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
  
  # === Save Settings if Required ===
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
  
  # === Prepare Output ===
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

# === Helper Function to Create Train-Test Indices ===

create_train_test_indices <- function(data, cv_type, k = 5, percentage = 20, number_folds = 10) {
  n <- nrow(data)
  indices <- 1:n
  if (cv_type == "k-fold") {
    folds <- sample(rep(1:k, length.out = n))
    train_list <- lapply(1:k, function(x) indices[folds != x])
    test_list <- lapply(1:k, function(x) indices[folds == x])
  } else if (cv_type == "loo") {
    train_list <- lapply(indices, function(x) indices[-x])
    test_list <- lapply(indices, function(x) x)
  } else if (cv_type == "lpo") {
    size <- floor((percentage / 100) * n)
    folds <- replicate(number_folds, sample(indices, size), simplify = FALSE)
    train_list <- lapply(folds, function(test) setdiff(indices, test))
    test_list <- folds
  } else {
    stop("Unsupported cv_type.")
  }
  
  return(list(train = train_list, test = test_list))
}

