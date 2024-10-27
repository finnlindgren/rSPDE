#' @noRd
set_prior <- function(prior, default_mean, default_precision = 0.1) {
  # Validate default parameters
  if (!is.numeric(default_mean) || length(default_mean) != 1) {
    stop("default_mean must be a single numeric value.")
  }
  if (!is.numeric(default_precision) || length(default_precision) != 1 || default_precision <= 0) {
    stop("default_precision must be a positive numeric scalar.")
  }

  # Return default prior if none is provided
  if (is.null(prior)) {
    return(list(mean = default_mean, precision = default_precision))
  }

  # Ensure prior only contains allowed elements
  allowed_elements <- c("mean", "precision")
  invalid_elements <- setdiff(names(prior), allowed_elements)
  if (length(invalid_elements) > 0) {
    warning(sprintf("Invalid elements in prior: %s. Only 'mean' and 'precision' are allowed.",
                    paste(invalid_elements, collapse = ", ")))
  }

  # Validate and set 'mean'
  if (!is.null(prior$mean)) {
    if (!is.numeric(prior$mean)) {
      stop("'mean' must be numeric.")
    }
    if (length(prior$mean) > 1) {
      warning("Vector provided for 'mean'; only the first element will be used.")
      prior$mean <- prior$mean[1]
    }
  } else {
    prior$mean <- default_mean  # Use default mean if not provided
  }

  # Validate and set 'precision'
  if (!is.null(prior$precision)) {
    if (!is.numeric(prior$precision) || any(prior$precision <= 0)) {
      stop("'precision' must contain positive numeric values.")
    }
    if (length(prior$precision) > 1) {
      warning("Vector provided for 'precision'; only the first element will be used.")
      prior$precision <- prior$precision[1]
    }
  } else {
    prior$precision <- default_precision  # Use default precision if not provided
  }
  return(prior)
}

#' Space-Time Random Fields via SPDE Approximation
#' 
#' `rspde.spacetime` computes a Finite Element Method (FEM) approximation of a 
#' Gaussian random field defined as the solution to the stochastic partial 
#' differential equation (SPDE):
#' \deqn{d u + \gamma(\kappa^2 + \rho\cdot\nabla - \Delta)^\alpha u = \sigma dW_C}
#' where \( C \) is a Whittle-Matérn covariance operator with smoothness parameter 
#' \(\beta\) and range parameter \(\kappa\). This function is designed to handle 
#' space-time random fields using either 1D spatial models or higher-dimensional 
#' FEM-based approaches.
#' 
#' @param mesh_space Spatial mesh for the FEM approximation, or a `metric_graph` 
#' object for handling models on metric graphs.
#' @param mesh_time Temporal mesh for the FEM approximation.
#' @param space_loc A vector of spatial locations for mesh nodes in 1D spatial models. 
#' This should be provided when `mesh_space` is not specified.
#' @param time_loc A vector of temporal locations for mesh nodes. This should be 
#' provided when `mesh_time` is not specified.
#' @param drift Logical value indicating whether the drift term should be included. 
#' If `FALSE`, the drift coefficient `\rho` is set to zero.
#' @param alpha Integer smoothness parameter \(\alpha\).
#' @param beta Integer smoothness parameter \(\beta\).
#' @param prior.kappa A list specifying the prior for the range parameter \(\kappa\). 
#' This list may contain two elements: `mean` and/or `precision`, both of which must 
#' be numeric scalars (numeric vectors of length 1). If `NULL`, default values will be used.
#' @param prior.sigma A list specifying the prior for the variance parameter \(\sigma\). 
#' This list may contain two elements: `mean` and/or `precision`, both of which must 
#' be numeric scalars. If `NULL`, default values will be used.
#' @param prior.rho A list specifying the prior for the drift coefficient \(\rho\). 
#' This list may contain two elements: `mean` and/or `precision`, both of which must 
#' be numeric scalars. If `NULL`, default values will be used. Will not be used if `drift = FALSE`.
#' @param prior.gamma A list specifying the prior for the weight \(\gamma\) in the SPDE 
#' operator. This list may contain two elements: `mean` and/or `precision`, both of which 
#' must be numeric scalars. If `NULL`, default values will be used.
#' @param shared_lib String specifying which shared library to use for the Cgeneric 
#' implementation. Options are "detect", "INLA", or "rSPDE". You may also specify the 
#' direct path to a .so (or .dll) file.
#' @param ... Additional arguments passed internally for configuration purposes.
#' @return An object of class `inla_rspde_spacetime` representing the FEM approximation of 
#' the space-time Gaussian random field.
#' @export
#' 
#' @examples
#' library(INLA)
#' library(MetricGraph)
#' graph <- metric_graph$new()
#' graph$build_mesh(h = 0.1)
#' graph$compute_fem()
#' 
#' # Define the time locations
#' time_loc <- seq(from = 0, to = 10, length.out = 11)
#' 
#' # Create the model
#' model <- rspde.spacetime(graph = graph, 
#'                          time_loc = time_loc, 
#'                          include_drift = FALSE, 
#'                          alpha = 2, 
#'                          beta = 1)
#' 

rspde.spacetime <- function(mesh_space = NULL,
                            mesh_time = NULL,
                            space_loc = NULL,
                            time_loc = NULL,
                            drift = TRUE,
                            alpha,
                            beta,
                            prior.kappa = NULL,
                            prior.sigma = NULL,
                            prior.rho = NULL,
                            prior.gamma = NULL,
                            shared_lib = "detect",
                            ...) {
  if (!is.null(mesh_space) && !is.null(space_loc)) {
    stop("Provide only one of 'mesh_space' or 'space_loc', not both.")
  }

  if (is.null(mesh_time) && is.null(time_loc)) {
    stop("You must provide either 'mesh_time' or 'time_loc'. Both cannot be NULL.")
  }

  if (is.null(alpha) || alpha %% 1 != 0) {
    stop("'alpha' must be provided and must be an integer.")
  }

  if (is.null(beta) || beta %% 1 != 0) {
    stop("'beta' must be provided and must be an integer.")
  }

  graph <- if (inherits(mesh_space, "metric_graph")) mesh_space else NULL
  if (!is.null(graph)) {
    mesh_space <- NULL
  }

  op <- spacetime.operators(
    mesh_space = mesh_space,
    mesh_time = mesh_time,
    space_loc = space_loc,
    time_loc = time_loc,
    graph = graph,
    kappa = prior.kappa$mean,
    sigma = prior.sigma$mean,
    gamma = prior.gamma$mean,
    rho = if (drift) prior.rho$mean else 0,
    alpha = alpha,
    beta = beta
  )

  prior.kappa <- set_prior(prior.kappa, op$kappa)
  prior.sigma <- set_prior(prior.sigma, op$sigma)
  prior.gamma <- set_prior(prior.gamma, op$gamma)
  prior.rho <- set_prior(prior.rho, op$rho)

  rspde_lib <- get_shared_library(shared_lib)

  list_args <- c(
    list(
      model = "inla_cgeneric_rspde_spacetime",
      shlib = rspde_lib,
      n = nrow(op$Q),
      d = op$d,
      Q = op$Q,
      n_Gtlist = length(op$Gtlist),
      n_Ctlist = length(op$Ctlist),
      n_B0list = length(op$B0list),
      n_M2list = length(op$M2list),
      n_M2list2 = length(op$M2list2),
      prior.kappa.mean = prior.kappa$mean,
      prior.kappa.precision = prior.kappa$precision,
      prior.sigma.mean = prior.sigma$mean,
      prior.sigma.precision = prior.sigma$precision,
      prior.gamma.mean = prior.gamma$mean,
      prior.gamma.precision = prior.gamma$precision,
      prior.rho.mean = prior.rho$mean,
      prior.rho.precision = prior.rho$precision,
      beta = beta,
      alpha = alpha,
      drift = as.integer(drift)
    ),
    # Single-level lists
    setNames(op$Gtlist, paste0("Gtlist", seq_along(op$Gtlist))),
    setNames(op$Ctlist, paste0("Ctlist", seq_along(op$Ctlist))),
    setNames(op$B0list, paste0("B0list", seq_along(op$B0list))),
    # Two-level M2list with lengths
    unlist(lapply(seq_along(op$M2list), function(i) {
      c(setNames(list(length(op$M2list[[i]])), paste0("n_M2list_", i)),
        setNames(op$M2list[[i]], paste0("M2list", i, seq_along(op$M2list[[i]]))))
    }), recursive = FALSE),
    # Two-level M2list2 with lengths, only if it has elements
    if (length(op$M2list2) > 0) unlist(lapply(seq_along(op$M2list2), function(i) {
      c(setNames(list(length(op$M2list2[[i]])), paste0("n_M2list2_", i)),
        setNames(op$M2list2[[i]], paste0("M2list2", i, seq_along(op$M2list2[[i]]))))
    }), recursive = FALSE)
  )

  model <- do.call(INLA::inla.cgeneric.define, list_args)

  model$A <- op$make_A
  model$drift <- drift
  model$prior.kappa <- prior.kappa
  model$prior.sigma <- prior.sigma
  model$prior.gamma <- prior.gamma
  model$prior.rho <- prior.rho
  model$d <- op$d
  
  class(model) <- c("inla_rspde_spacetime", class(model))

  return(model)
}


# make.L <- function(n,kappa,Glist){
#     if(n > length(Glist)){
#         stop("Glist too short.")
#     }
#     L <- 0
#     for(k in 0:n){
#         L <- L + choose(n,k)*kappa^(2*(n-k))*Glist[[k+1]]
#     }
#     return(L)
# }

# # Initialize Q with the first term
# Q <- make.L(beta, kappa, Gtlist) + 2 * gamma * make.L(beta + alpha, kappa, B0list)

# # Add contributions from alpha for Q
# for (k in 0:alpha) {
#     # Add the gamma^2 term
#     Q <- Q + gamma^2 * choose(alpha, k) * rho^(2 * k) * make.L(beta + 2 * (alpha - k), kappa, Ctlist[(k + 1):length(Ctlist)])

#     # Add M2 list terms for 1D and 2D cases
#     if (d == 2) {
#         M2x <- make.L(beta + alpha - k, kappa, M2list[[k + 1]])
#         M2y <- make.L(beta + alpha - k, kappa, M2list2[[k + 1]])
#         Q <- Q - 0.5 * gamma * choose(alpha, k) * (1 - (-1)^k) * rho[1]^(k) * (M2x + t(M2x))
#         Q <- Q - 0.5 * gamma * choose(alpha, k) * (1 - (-1)^k) * rho[2]^(k) * (M2y + t(M2y))
#     } else {
#         M2 <- make.L(beta + alpha - k, kappa, M2list[[k + 1]])
#         Q <- Q - 0.5 * (-1)^(floor(k / 2)) * gamma * choose(alpha, k) * (1 - (-1)^k) * rho^(k) * (M2 + t(M2))
#     }
# }