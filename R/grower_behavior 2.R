#' Delineate management units from a host raster
#'
#' Groups raster cells into management units (farms/fields) for which a single
#' grower makes a joint treatment decision. Three methods are supported:
#' 'raster' (contiguous host patches via terra::patches), 'polygon' (rasterise
#' supplied field/farm boundaries), and 'grid' (regular square grid).
#'
#' @param host_raster SpatRaster with host density values (NA or 0 = non-host)
#' @param method character; one of 'raster', 'polygon', 'grid'
#' @param polygon_file character; path to a vector file (shapefile or GeoPackage)
#'   with field/farm boundaries. Required when method = 'polygon'.
#' @param grid_size numeric; side length of square management units in map units
#'   (e.g. metres). Required when method = 'grid'.
#' @param min_cells integer; minimum number of host cells a patch must contain
#'   to be retained as a management unit (patches smaller than this are merged
#'   into their nearest neighbor). Default 1.
#'
#' @importFrom terra patches rasterize ext res values as.matrix rowFromCell
#'   colFromCell xyFromCell nrow ncol app rast vect
#' @return A named list with:
#'   \item{management_unit_matrix}{integer matrix of unit IDs (0 = non-unit)}
#'   \item{management_unit_table}{data.frame with unit_id, n_cells,
#'     centroid_x, centroid_y}
#' @export
delineate_management_units <- function(host_raster,
                                        method      = "raster",
                                        polygon_file = NULL,
                                        grid_size    = NULL,
                                        min_cells    = 1) {

  if (!method %in% management_unit_method_list) {
    stop(management_unit_method_error)
  }
  if (method == "grid" && is.null(grid_size)) {
    stop(grid_size_missing_error)
  }
  if (method == "polygon" && is.null(polygon_file)) {
    stop(polygon_file_missing_error)
  }

  # Binary host mask (1 where host present, NA elsewhere)
  host_mask <- host_raster[[1]]
  host_mask[host_mask <= 0] <- NA

  if (method == "raster") {
    # Label contiguous patches of host cells
    patches_rast <- terra::patches(host_mask, directions = 8, zeroAsNA = TRUE)

  } else if (method == "polygon") {
    polys <- terra::vect(polygon_file)
    # Rasterise polygon IDs onto the host raster grid
    patches_rast <- terra::rasterize(polys, host_raster[[1]], field = 1:nrow(polys))
    # Mask to host-only cells
    patches_rast <- terra::mask(patches_rast, host_mask)
    # Re-label to contiguous integers
    uid  <- sort(unique(terra::values(patches_rast, na.rm = TRUE)))
    vals <- terra::values(patches_rast)
    new_vals <- match(vals, uid)   # NA stays NA
    terra::values(patches_rast) <- new_vals

  } else {  # "grid"
    r <- terra::res(host_raster[[1]])
    # Number of host cells per grid square
    fact_x <- max(1L, round(grid_size / r[1]))
    fact_y <- max(1L, round(grid_size / r[2]))
    # Create grid IDs by row/col block
    nr <- terra::nrow(host_raster[[1]])
    nc <- terra::ncol(host_raster[[1]])
    row_ids <- ceiling(seq_len(nr) / fact_y)
    col_ids <- ceiling(seq_len(nc) / fact_x)
    id_mat  <- outer(row_ids, col_ids, function(r, c) (r - 1) * max(col_ids) + c)
    patches_rast <- host_raster[[1]]
    # terra::values() fills in row-major order; as.integer() on a matrix is
    # column-major, so transpose first to get the correct spatial layout.
    terra::values(patches_rast) <- as.integer(t(id_mat))
    patches_rast <- terra::mask(patches_rast, host_mask)
  }

  # Replace NAs with 0 (background)
  patches_vals <- terra::values(patches_rast)
  patches_vals[is.na(patches_vals)] <- 0L
  terra::values(patches_rast) <- as.integer(patches_vals)

  unit_matrix <- terra::as.matrix(patches_rast, wide = TRUE)
  unit_matrix[is.na(unit_matrix)] <- 0L
  storage.mode(unit_matrix) <- "integer"

  # Enforce min_cells: relabel small patches to 0
  unit_ids <- sort(unique(as.integer(unit_matrix[unit_matrix > 0])))
  counts   <- tabulate(unit_matrix[unit_matrix > 0])
  # tabulate uses indices starting at 1, so align with unit_ids
  count_vec <- vapply(unit_ids, function(id) sum(unit_matrix == id), integer(1))
  small_ids <- unit_ids[count_vec < min_cells]
  if (length(small_ids) > 0) {
    unit_matrix[unit_matrix %in% small_ids] <- 0L
  }

  # Rebuild contiguous integer IDs after removal
  remaining_ids <- sort(unique(as.integer(unit_matrix[unit_matrix > 0])))
  id_map <- setNames(seq_along(remaining_ids), remaining_ids)
  new_mat <- unit_matrix
  for (old_id in remaining_ids) {
    new_mat[unit_matrix == old_id] <- id_map[as.character(old_id)]
  }
  unit_matrix <- new_mat

  # Build lookup table with centroids.
  # xyFromCell() numbers cells in row-major order (cell 1 = row 1 col 1,
  # cell 2 = row 1 col 2, …).  as.integer(matrix) is column-major, so we
  # must NOT use a flat vector index to look up xy rows.  Instead we use
  # which(arr.ind = TRUE) to get explicit [row, col] positions and convert
  # them to terra cell numbers before calling xyFromCell.
  unit_ids_final <- sort(unique(as.integer(unit_matrix[unit_matrix > 0])))
  nc_rast <- terra::ncol(host_raster[[1]])

  table_rows <- lapply(unit_ids_final, function(uid) {
    rc    <- which(unit_matrix == uid, arr.ind = TRUE)   # [row, col] pairs
    cells <- (rc[, 1] - 1L) * nc_rast + rc[, 2]         # terra row-major cell numbers
    xy_uid <- terra::xyFromCell(host_raster[[1]], cells)
    data.frame(
      unit_id    = uid,
      n_cells    = nrow(rc),
      centroid_x = mean(xy_uid[, 1]),
      centroid_y = mean(xy_uid[, 2])
    )
  })
  unit_table <- do.call(rbind, table_rows)

  return(list(
    management_unit_matrix = unit_matrix,
    management_unit_table  = unit_table
  ))
}

# ──────────────────────────────────────────────────────────────────────────────

#' Assign grower behavioral types to management units
#'
#' Allocates each management unit to one of K grower behavioral types using
#' one of three spatial structures: independent random assignment, spatially
#' correlated (clustered) assignment via a Gaussian random field, or
#' empirically constrained assignment based on an adoption-probability raster.
#'
#' @param management_units integer matrix of management unit IDs (0 = background)
#'   as returned by \code{delineate_management_units()$management_unit_matrix}
#' @param type_names character vector of grower type names, e.g.
#'   \code{c("early_adopter", "responsive", "non_adopter")}
#' @param type_probs numeric vector of prior probabilities for each type.
#'   Must sum to 1 and be the same length as \code{type_names}.
#' @param spatial_structure one of \code{"random"}, \code{"clustered"},
#'   \code{"empirical"}
#' @param cluster_range numeric; correlation range parameter in map units
#'   (metres) for the Gaussian random field. Required for
#'   \code{spatial_structure = "clustered"}.
#' @param cluster_model one of \code{"gaussian"} or \code{"exponential"};
#'   the covariance model used for the latent field.
#' @param empirical_raster SpatRaster; each cell gives the probability of
#'   being the *first* type in \code{type_names}. Required for
#'   \code{spatial_structure = "empirical"}.
#' @param host_raster SpatRaster used only to obtain spatial metadata
#'   (extent, resolution, CRS) when \code{spatial_structure = "clustered"}.
#'   Ignored for other structures.
#' @param random_seed integer or NULL; passed to \code{set.seed()} for
#'   reproducibility.
#'
#' @return integer matrix (same dimensions as \code{management_units}) with
#'   values \code{1..K} indicating the grower type index. Background cells
#'   (unit ID 0) are set to 0.
#' @export
assign_grower_types <- function(management_units,
                                 type_names,
                                 type_probs,
                                 spatial_structure = "random",
                                 cluster_range     = NULL,
                                 cluster_model     = "gaussian",
                                 empirical_raster  = NULL,
                                 host_raster       = NULL,
                                 random_seed       = NULL) {

  # --- Input validation -------------------------------------------------------
  check <- behavior_type_checks(type_names, type_probs, spatial_structure,
                                  cluster_model)
  if (!check$checks_passed) stop(check$failed_check)

  if (!is.null(random_seed)) set.seed(random_seed)

  K         <- length(type_names)
  unit_ids  <- sort(unique(as.integer(management_units[management_units > 0])))
  N         <- length(unit_ids)

  # --- Assign types -----------------------------------------------------------
  if (spatial_structure == "random") {
    type_assignments <- sample.int(K, size = N, replace = TRUE, prob = type_probs)

  } else if (spatial_structure == "clustered") {
    if (is.null(cluster_range) || cluster_range <= 0) {
      stop("cluster_range must be a positive number when spatial_structure = 'clustered'")
    }
    # Simulate a latent Gaussian random field at unit centroids.
    # We use a simple exponential or Gaussian semivariogram applied via
    # distance-decay weighting — no external spatial-stats package required.
    # For each unit, compute a weighted average of K independent N(0,1) fields
    # using a covariance kernel; then threshold to type boundaries.

    # Get centroid coordinates for each unit
    nr <- nrow(management_units)
    nc <- ncol(management_units)
    # Approximate centroids from matrix position (row/col center)
    centroids <- vapply(unit_ids, function(uid) {
      idx <- which(management_units == uid, arr.ind = TRUE)
      c(mean(idx[, 2]), mean(nr - idx[, 1] + 1))  # (col, row) as (x, y)
    }, numeric(2))
    centroids <- t(centroids)  # N x 2

    # Pairwise distances (in cell units, scale by notional resolution = 1)
    dist_mat <- as.matrix(dist(centroids))

    # Covariance kernel
    if (cluster_model == "gaussian") {
      cov_mat <- exp(-(dist_mat / cluster_range)^2)
    } else {  # exponential
      cov_mat <- exp(-dist_mat / cluster_range)
    }
    # Regularise for Cholesky
    cov_mat <- cov_mat + diag(1e-6, N)

    # Cholesky factorisation for correlated draws
    L <- tryCatch(chol(cov_mat), error = function(e) {
      # Fallback: nearest positive-definite via eigendecomposition
      eig <- eigen(cov_mat, symmetric = TRUE)
      eig$values[eig$values < 1e-8] <- 1e-8
      L_approx <- eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors)
      chol(L_approx + diag(1e-6, N))
    })

    # Draw one correlated standard-normal field
    z <- as.numeric(t(L) %*% rnorm(N))

    # Convert to type probabilities via cumulative normal thresholds
    cum_probs <- c(0, cumsum(type_probs))
    thresholds <- stats::qnorm(cum_probs)  # -Inf, q1, q2, ..., Inf
    type_assignments <- as.integer(cut(z, breaks = thresholds, labels = seq_len(K)))

  } else {  # "empirical"
    if (is.null(empirical_raster)) {
      stop("empirical_raster must be supplied when spatial_structure = 'empirical'")
    }
    # Extract per-unit mean of empirical_raster (probability of type 1)
    nr <- nrow(management_units)
    nc <- ncol(management_units)
    emp_mat <- as.matrix(empirical_raster[[1]], wide = TRUE)
    if (!identical(dim(emp_mat), dim(management_units))) {
      stop("empirical_raster must have the same dimensions as management_units")
    }

    type_assignments <- vapply(unit_ids, function(uid) {
      idx <- which(management_units == uid)
      p1  <- mean(emp_mat[idx], na.rm = TRUE)
      p1  <- max(0, min(1, p1))
      # Remaining probability distributed proportionally among other types
      if (K == 1) return(1L)
      other_sum   <- sum(type_probs[-1])
      other_probs <- if (other_sum > 0) type_probs[-1] / other_sum * (1 - p1) else
                       rep((1 - p1) / (K - 1), K - 1)
      p_vec <- c(p1, other_probs)
      as.integer(sample.int(K, size = 1, prob = p_vec))
    }, integer(1))
  }

  # --- Build output matrix ---------------------------------------------------
  out_mat <- matrix(0L, nrow = nrow(management_units), ncol = ncol(management_units))
  for (i in seq_along(unit_ids)) {
    out_mat[management_units == unit_ids[i]] <- type_assignments[i]
  }
  storage.mode(out_mat) <- "integer"
  return(out_mat)
}

# ──────────────────────────────────────────────────────────────────────────────

#' Grower perception and treatment decision function
#'
#' At a management decision point, reads the current infection state, applies a
#' per-type perception filter (detection probability), and returns a treatment
#' decision matrix that can be passed directly to the PoPS treatment system.
#'
#' Prevalence is computed in two passes. Pass 1 aggregates infected and host
#' cell counts for every management unit (field). Pass 2 applies the decision
#' model for each unit, optionally expanding the observation neighborhood to
#' include all units whose centroid falls within \code{spatial_radius} map units
#' of the focal unit's centroid before applying detection and the decision
#' threshold. This follows Chris Jones's recommendation to track infected area
#' per grower ID as columns in the management unit table.
#'
#' @param infected_matrix integer matrix; current infection counts per cell
#' @param host_matrix integer matrix; total host counts per cell
#' @param grower_type_matrix integer matrix; grower type index per cell
#'   (0 = background), as returned by \code{assign_grower_types()}
#' @param management_unit_matrix integer matrix; unit IDs per cell (0 = background),
#'   as returned by \code{delineate_management_units()$management_unit_matrix}
#' @param management_unit_table data.frame; lookup table with columns
#'   \code{unit_id}, \code{n_cells}, \code{centroid_x}, \code{centroid_y}
#' @param behavior_params named list; one entry per type name (matching the
#'   integer index in \code{grower_type_matrix}), each a list with:
#'   \describe{
#'     \item{detection_prob}{numeric in [0,1]; probability that an infected cell
#'       within the observation neighborhood is detected}
#'     \item{willingness_to_treat}{numeric in [0,1]; probability of deciding
#'       to treat given perceived prevalence exceeds the threshold}
#'     \item{treatment_efficacy}{numeric in [0,1]; fractional reduction in
#'       infection applied by the treatment}
#'     \item{decision_threshold}{numeric in [0,1]; perceived prevalence above
#'       which the grower considers treating}
#'     \item{spatial_radius}{numeric (optional); radius in map units (same CRS
#'       as centroid coordinates) within which the grower observes infection.
#'       Growers aggregate infected and host counts across all units whose
#'       centroid falls within this distance. If \code{NULL} or absent, only
#'       the grower's own unit is observed.}
#'   }
#' @param output character; \code{"binary"} returns a 0/1 decision matrix;
#'   \code{"continuous"} returns treatment efficacy values per cell
#' @param random_seed integer or NULL
#'
#' @return numeric matrix with the same dimensions as \code{host_matrix};
#'   0 = no treatment, 1 (or efficacy) = treat. Suitable as a PoPS
#'   \code{treatment_map}. The \code{management_unit_table} attribute on the
#'   return value includes \code{infected_cells}, \code{host_cells},
#'   \code{prevalence}, \code{perceived_prevalence}, and \code{treated} columns
#'   populated for this decision step.
#' @export
grower_decision <- function(infected_matrix,
                             host_matrix,
                             grower_type_matrix,
                             management_unit_matrix,
                             management_unit_table,
                             behavior_params,
                             output      = "binary",
                             random_seed = NULL) {

  if (!output %in% decision_output_list) stop(decision_output_error)
  if (!is.null(random_seed)) set.seed(random_seed)

  nr      <- nrow(host_matrix)
  nc      <- ncol(host_matrix)
  out     <- matrix(0.0, nrow = nr, ncol = nc)
  n_units <- nrow(management_unit_table)
  unit_ids <- management_unit_table$unit_id

  # ── Pass 1: Aggregate per-unit infection counts ───────────────────────────────
  # Infected area per grower ID stored as columns in management_unit_table
  # (Chris Jones feedback); used in Pass 2 for neighborhood aggregation.
  management_unit_table$infected_cells <- 0L
  management_unit_table$host_cells     <- 0L
  management_unit_table$prevalence     <- 0.0

  for (i in seq_len(n_units)) {
    uid        <- unit_ids[i]
    unit_cells <- which(management_unit_matrix == uid, arr.ind = TRUE)
    if (nrow(unit_cells) == 0) next

    unit_infected <- sum(infected_matrix[unit_cells], na.rm = TRUE)
    unit_hosts    <- sum(host_matrix[unit_cells],     na.rm = TRUE)

    management_unit_table$infected_cells[i] <- as.integer(unit_infected)
    management_unit_table$host_cells[i]     <- as.integer(unit_hosts)
    management_unit_table$prevalence[i]     <-
      if (unit_hosts > 0) unit_infected / unit_hosts else 0.0
  }

  # ── Pass 2: Decision loop with optional spatial_radius neighborhood ──────────
  management_unit_table$perceived_prevalence <- NA_real_
  management_unit_table$treated              <- FALSE

  for (i in seq_len(n_units)) {
    uid        <- unit_ids[i]
    unit_cells <- which(management_unit_matrix == uid, arr.ind = TRUE)
    if (nrow(unit_cells) == 0) next

    # Determine grower type for this unit (modal type across its cells)
    type_vals <- grower_type_matrix[unit_cells]
    type_val  <- as.integer(names(sort(table(type_vals[type_vals > 0]),
                                       decreasing = TRUE)[1]))
    if (length(type_val) == 0 || is.na(type_val)) next

    params <- behavior_params[[type_val]]
    if (is.null(params)) next

    # ── 1. Determine observation neighborhood ──────────────────────────────────
    spatial_radius <- params$spatial_radius %||% NULL

    if (!is.null(spatial_radius) && is.numeric(spatial_radius) &&
        spatial_radius > 0 && n_units > 1) {

      # Euclidean distance from this unit's centroid to every other centroid
      dx   <- management_unit_table$centroid_x - management_unit_table$centroid_x[i]
      dy   <- management_unit_table$centroid_y - management_unit_table$centroid_y[i]
      dist <- sqrt(dx^2 + dy^2)

      # Include units whose centroid is within spatial_radius (self always included)
      nbrs         <- which(dist <= spatial_radius)
      obs_infected <- sum(management_unit_table$infected_cells[nbrs], na.rm = TRUE)
      obs_host     <- sum(management_unit_table$host_cells[nbrs],     na.rm = TRUE)

    } else {
      # No radius: observe own unit only
      obs_infected <- management_unit_table$infected_cells[i]
      obs_host     <- management_unit_table$host_cells[i]
    }

    # ── 2. Apply detection probability (binomial thinning) ──────────────────────
    detection_prob <- params$detection_prob %||% 1.0
    if (obs_host > 0 && detection_prob < 1) {
      observed_infected <- stats::rbinom(1, size = obs_infected,
                                          prob = detection_prob)
    } else {
      observed_infected <- obs_infected * detection_prob
    }
    perceived_prevalence <- if (obs_host > 0) observed_infected / obs_host else 0
    management_unit_table$perceived_prevalence[i] <- perceived_prevalence

    # ── 3. Decision rule ────────────────────────────────────────────────────────
    decision_threshold   <- params$decision_threshold   %||% 0.0
    willingness_to_treat <- params$willingness_to_treat %||% 0.0
    treatment_efficacy   <- params$treatment_efficacy   %||% 1.0

    treats <- FALSE
    if (perceived_prevalence > decision_threshold) {
      treats <- stats::runif(1) < willingness_to_treat
    }
    management_unit_table$treated[i] <- treats

    if (treats) {
      efficacy_val <- if (output == "binary") 1.0 else treatment_efficacy
      out[unit_cells] <- efficacy_val
    }
  }

  attr(out, "management_unit_table") <- management_unit_table
  return(out)
}

# ──────────────────────────────────────────────────────────────────────────────

#' Parse and validate behavior module configuration
#'
#' Analogous to \code{configuration()} for the core model. Accepts either a
#' named list or a path to a YAML file and returns a validated
#' \code{behavior_config} list.
#'
#' The function handles the following YAML fields:
#' \describe{
#'   \item{use_behavior_module}{logical toggle (default FALSE); mirrors
#'     \code{use_movements} etc. in the core config. The legacy name
#'     \code{enable_behavior_module} is accepted as an alias.}
#'   \item{grower_id_file}{character; path to a GeoTIFF with integer grower/
#'     management-unit IDs (0 = background). Produced by
#'     \code{delineate_management_units()}.}
#'   \item{grower_params_file}{character; path to a CSV with a \code{grower_id}
#'     column plus parameter columns
#'     (\code{detection_prob}, \code{willingness_to_treat},
#'     \code{treatment_efficacy}, \code{decision_threshold}). When supplied, its
#'     values override \code{behavior_params} for matching grower IDs.}
#'   \item{behavior_decision_dates}{character vector of YYYY-MM-DD dates at
#'     which growers make treatment decisions.}
#'   \item{behavior_params}{list; one entry per grower type (used as global
#'     defaults / fallback when \code{grower_params_file} is absent or
#'     incomplete).}
#' }
#'
#' @param config_input named list or character path to a YAML file
#' @param start_date simulation start date (YYYY-MM-DD) for date range checks
#' @param end_date simulation end date (YYYY-MM-DD) for date range checks
#'
#' @return named list with validated behavior parameters; includes
#'   \code{$failure} if any check fails (mirrors \code{configuration()})
#' @importFrom yaml yaml.load_file
#' @export
behavior_configuration <- function(config_input, start_date, end_date) {
  if (is.character(config_input)) {
    bc <- yaml::yaml.load_file(config_input)
  } else {
    bc <- config_input
  }

  bc$failure <- NULL

  # ── Toggle ────────────────────────────────────────────────────────────────
  # Canonical field name follows PoPS convention: use_behavior_module.
  # Accept the legacy alias enable_behavior_module for backward compatibility.
  if (is.null(bc$use_behavior_module) && !is.null(bc$enable_behavior_module)) {
    bc$use_behavior_module <- bc$enable_behavior_module
  }
  bc$use_behavior_module <- bc$use_behavior_module %||% FALSE

  # ── Spatial data files ────────────────────────────────────────────────────
  bc$grower_id_file    <- bc$grower_id_file    %||% ""
  bc$grower_params_file <- bc$grower_params_file %||% ""

  # ── Decision dates ────────────────────────────────────────────────────────
  # Merge behavior_decision_dates (new name) and legacy decision_dates.
  if (is.null(bc$behavior_decision_dates) && !is.null(bc$decision_dates)) {
    bc$behavior_decision_dates <- bc$decision_dates
  }
  if (!is.null(bc$behavior_decision_dates)) {
    date_check <- decision_date_checks(bc$behavior_decision_dates, start_date, end_date)
    if (!date_check$checks_passed) {
      bc$failure <- date_check$failed_check
      return(bc)
    }
  }

  # ── Type assignment parameters ────────────────────────────────────────────
  # type_names / type_probs / spatial_structure are only required when
  # assign_grower_types() will be called (pre-processing step in R); they are
  # not needed by the C++ behavior module itself, so only validate if present.
  if (!is.null(bc$type_names)) {
    type_check <- behavior_type_checks(
      bc$type_names, bc$type_probs,
      bc$spatial_structure %||% "random",
      bc$cluster_model %||% "gaussian"
    )
    if (!type_check$checks_passed) {
      bc$failure <- type_check$failed_check
      return(bc)
    }
  }

  # ── Per-type behavior parameters (global defaults) ───────────────────────
  # Set defaults for any missing sub-params in each entry.
  if (!is.null(bc$behavior_params)) {
    for (i in seq_along(bc$behavior_params)) {
      p <- bc$behavior_params[[i]]
      p$detection_prob       <- p$detection_prob       %||% 1.0
      p$willingness_to_treat <- p$willingness_to_treat %||% 0.5
      p$treatment_efficacy   <- p$treatment_efficacy   %||% 1.0
      p$decision_threshold   <- p$decision_threshold   %||% 0.0
      bc$behavior_params[[i]] <- p
    }
  }

  return(bc)
}

# Null-coalescing operator (avoids dependency on rlang)
`%||%` <- function(a, b) if (!is.null(a)) a else b
