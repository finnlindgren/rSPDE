#' @noRd
#' @noRd
set_prior <- function(prior, default_mean, default_precision, p = 1) {
  # Validate default parameters
  if (!is.numeric(default_mean) || length(default_mean) != p) {
    stop(paste("default_mean must be a numeric vector of length equal to",p,"."))
  }
  if (!is.numeric(default_precision) || length(default_precision) != p || any(default_precision <= 0)) {
    stop(paste("default_precision must be a positive numeric vector of length equal to",p,"."))
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
    if (!is.numeric(prior$mean) || length(prior$mean) != p) {
      stop(sprintf("'mean' must be a numeric vector of length %d.", p))
    }
  } else {
    prior$mean <- default_mean  # Use default mean if not provided
  }

  # Validate and set 'precision'
  if (!is.null(prior$precision)) {
    if (!is.numeric(prior$precision) || length(prior$precision) != p || any(prior$precision <= 0)) {
      stop(sprintf("'precision' must be a positive numeric vector of length %d.", p))
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
#' where \eqn{C} is a Whittle-Matérn covariance operator with smoothness parameter
#' \eqn{\beta} and range parameter \eqn{\kappa}. This function is designed to handle
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
#' If `FALSE`, the drift coefficient \eqn{\rho} is set to zero.
#' @param alpha Integer smoothness parameter \eqn{\alpha}.
#' @param beta Integer smoothness parameter \eqn{\beta}.
#' @param prior.kappa A list specifying the prior for the range parameter \eqn{\kappa}.
#' This list may contain two elements: `mean` and/or `precision`, both of which must
#' be numeric scalars (numeric vectors of length 1). The precision refers to the prior
#' on \eqn{\log(\kappa)}. If `NULL`, default values will be used.
#' @param prior.sigma A list specifying the prior for the variance parameter \eqn{\sigma}.
#' This list may contain two elements: `mean` and/or `precision`, both of which must
#' be numeric scalars. The precision refers to the prior on \eqn{\log(\sigma)}. If `NULL`,
#' default values will be used.
#' @param prior.rho A list specifying the prior for the drift coefficient \eqn{\rho}.
#' This list may contain two elements: `mean` and/or `precision`, both of which must
#' be numeric scalars if dimension is one, and numeric vectors of length 2 if dimension is 2.
#' The precision applies directly to \eqn{\rho} without log transformation.
#' If `NULL`, default values will be used. Will not be used if `drift = FALSE`.
#' @param prior.gamma A list specifying the prior for the weight \eqn{\gamma} in the SPDE
#' operator. This list may contain two elements: `mean` and/or `precision`, both of which
#' must be numeric scalars. The precision refers to the prior on \eqn{\log(\gamma)}. If `NULL`,
#' default values will be used.
#' @param prior.precision A precision matrix for \eqn{\log(\kappa), \log(\sigma), \log(\gamma), \rho}. This matrix replaces the precision
#' element from `prior.kappa`, `prior.sigma`, `prior.gamma`, and `prior.rho` respectively. For dimension 1 `prior.precision` must be a 4x4 matrix. For dimension 2, \eqn{\rho} is a vector of length 2, so in this case `prior.precision` must be a 5x5 matrix. If `NULL`, a diagonal precision matrix with default values will be used.
#' @param shared_lib String specifying which shared library to use for the Cgeneric
#' implementation. Options are "detect", "INLA", or "rSPDE". You may also specify the
#' direct path to a .so (or .dll) file.
#' @param debug Logical value indicating whether to enable INLA debug mode.
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
#' model <- rspde.spacetime(mesh_space = graph,
#'                          time_loc = time_loc,
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
                            prior.precision = NULL,
                            shared_lib = "detect",
                            debug = FALSE,
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
    rho = prior.rho$mean,
    alpha = alpha,
    beta = beta
  )

  default_precision <- 0.1
  default_precision_rho <- if(op$d == 1) 0.1 else c(0.1, 0.1)

  prior.kappa <- set_prior(prior.kappa, op$kappa, default_precision, p = 1)
  prior.sigma <- set_prior(prior.sigma, op$sigma, default_precision, p = 1)
  prior.gamma <- set_prior(prior.gamma, op$gamma, default_precision, p = 1)
  prior.rho <- set_prior(prior.rho, op$rho, default_precision_rho, p = op$d)

  # Set default precision matrix if prior.precision is NULL
  if (is.null(prior.precision)) {
    if (drift) {
      # Include prior.rho$precision if drift is TRUE
      prior.precision <- diag(c(
        prior.kappa$precision,
        prior.sigma$precision,
        prior.gamma$precision,
        prior.rho$precision
      ))
    } else {
      # Exclude prior.rho$precision if drift is FALSE
      prior.precision <- diag(c(
        prior.kappa$precision,
        prior.sigma$precision,
        prior.gamma$precision
      ))
    }
  } else {
    # Check matrix dimensions based on op$d and drift
    expected_dim <- if (op$d == 1) {
      if (drift) c(4, 4) else c(3, 3)
    } else if (op$d == 2) {
      if (drift) c(5, 5) else c(3, 3)
    } else stop("dimension must be 1 or 2.")

    if (!is.matrix(prior.precision) || !all(dim(prior.precision) == expected_dim)) {
      stop(sprintf("prior.precision must be a %dx%d matrix.", expected_dim[1], expected_dim[2]))
    }
  }

  rspde_lib <- get_shared_library(shared_lib)

  Cmatrix  <-  as(as(as(op$Q, "dMatrix"), "generalMatrix"), "TsparseMatrix")
  ii <- Cmatrix@i
  Cmatrix@i <- Cmatrix@j
  Cmatrix@j <- ii
  idx <- which(Cmatrix@i <= Cmatrix@j)
  Cmatrix@i <- Cmatrix@i[idx]
  Cmatrix@j <- Cmatrix@j[idx]
  Cmatrix@x <- Cmatrix@x[idx]

  list_args <- c(
    list(
      model = "inla_cgeneric_rspde_spacetime_model",
      shlib = rspde_lib,
      n = nrow(op$Q),
      debug = debug,
      d = as.integer(op$d),
      Q = Cmatrix,
      n_Gtlist = length(op$Gtlist),
      n_Ctlist = length(op$Ctlist),
      n_B0list = length(op$B0list),
      n_M2list = length(op$M2list),
      n_M2list2 = length(op$M2list2),
      prior.kappa.mean = prior.kappa$mean,
      prior.sigma.mean = prior.sigma$mean,
      prior.gamma.mean = prior.gamma$mean,
      prior.rho.mean = prior.rho$mean,
      beta = beta,
      alpha = alpha,
      drift = as.integer(drift),
      prior.precision = prior.precision
    ),

    # Single-level lists, each element added separately
    setNames(op$Gtlist, paste0("Gtlist", seq_along(op$Gtlist))),
    setNames(op$Ctlist, paste0("Ctlist", seq_along(op$Ctlist))),
    setNames(op$B0list, paste0("B0list", seq_along(op$B0list))),

    # Flattened two-level M2list
    unlist(lapply(seq_along(op$M2list), function(i) {
      c(
        # Add the length of each first-level M2list element
        setNames(list(length(op$M2list[[i]])), paste0("n_M2list_", i)),

        # Add each second-level element with a unique name
        setNames(op$M2list[[i]], paste0("M2list", i, "_", seq_along(op$M2list[[i]])))
      )
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
  model$prior.precision <- prior.precision
  if(!is.null(op$mesh_space)){
    model$mesh <- op$mesh_space
  } else{
    model$mesh <- graph
  }
  model$time_mesh <- op$mesh_time

  class(model) <- c("inla_rspde_spacetime", class(model))

  return(model)
}



#' @title rSPDE space time inlabru mapper
#' @name bru_get_mapper.inla_rspde_spacetime
#' @param model An `inla_rspde_spacetime` object for which to construct or extract a mapper
#' @param \dots Arguments passed on to other methods
#' @rdname bru_get_mapper.inla_rspde_spacetime
#' @rawNamespace if (getRversion() >= "3.6.0") {
#'   S3method(inlabru::bru_get_mapper, inla_rspde_spacetime)
#'   S3method(inlabru::bru_mapper, metric_graph)
#'   S3method(inlabru::ibm_n, bru_mapper_metric_graph)
#'   S3method(inlabru::ibm_values, bru_mapper_metric_graph)
#'   S3method(inlabru::ibm_jacobian, bru_mapper_metric_graph)
#' }
#'

bru_get_mapper.inla_rspde_spacetime <- function(model, ...) {
  stopifnot(requireNamespace("inlabru"))
  inlabru::bru_mapper_multi(list(
    space = if(inherits(model[["mesh"]], c("fm_mesh_1d", "inla.mesh.1d"))){
      inlabru::bru_mapper(model[["mesh"]], indexed = TRUE)
    } else{
      inlabru::bru_mapper(model[["mesh"]])
    },
    time = inlabru::bru_mapper(model[["time_mesh"]], indexed = TRUE)
  ))
}

#' @noRd
bru_mapper.metric_graph <- function(mesh, ...) {
  mapper <- list(mesh = mesh)
  inlabru::bru_mapper_define(mapper, new_class = "bru_mapper_metric_graph")
}

#' @noRd
ibm_n.bru_mapper_metric_graph <- function(mapper, ...) {
  nrow(mapper[["mesh"]][["mesh"]][["VtE"]])
}

#' @noRd
ibm_values.bru_mapper_metric_graph <- function(mapper, ...) {
  seq_len(nrow(mapper[["mesh"]][["mesh"]][["VtE"]]))
}

#' @noRd
ibm_jacobian.bru_mapper_metric_graph <- function(mapper, input, ...) {
  if (is.null(input)) {
    return(Matrix::Matrix(0, 0, inlabru::ibm_n(mapper)))
  }
  mapper[["mesh"]][["fem_basis"]](input)
}
