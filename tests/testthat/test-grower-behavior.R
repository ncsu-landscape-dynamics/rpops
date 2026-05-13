## Tests for the grower/manager behavior module
##
## Covers:
##   delineate_management_units()
##   assign_grower_types()
##   grower_decision()
##   behavior_configuration()
##
## Functions removed in the architectural pivot to C++-native behavior
## (pops_with_behavior, extract_host_pool_state, build_treatment_map,
## update_grower_trust) no longer have test coverage here.  Tests for the
## C++ BehaviorModule will live in test-pops.r once behavior.hpp is wired
## into pops.cpp and pops_model_cpp().

library(testthat)
library(terra)

# -- Shared test fixtures -----------------------------------------------------

make_host_raster <- function() {
  r <- terra::rast(nrows = 4, ncols = 4, xmin = 0, xmax = 4, ymin = 0, ymax = 4)
  terra::values(r) <- c(
    1, 1, 0, 0,
    1, 1, 0, 0,
    0, 0, 1, 1,
    0, 0, 1, 1
  )
  r
}

make_decision_setup <- function() {
  nr <- 4; nc <- 4
  host_mat     <- matrix(c(1, 1, 0, 0,
                           1, 1, 0, 0,
                           0, 0, 1, 1,
                           0, 0, 1, 1), nrow = nr, byrow = TRUE)
  infected_mat <- matrix(c(1, 0, 0, 0,
                           0, 0, 0, 0,
                           0, 0, 1, 0,
                           0, 0, 0, 0), nrow = nr, byrow = TRUE)
  # Two management units: top-left block = unit 1, bottom-right = unit 2
  unit_mat <- matrix(c(1, 1, 0, 0,
                       1, 1, 0, 0,
                       0, 0, 2, 2,
                       0, 0, 2, 2), nrow = nr, byrow = TRUE)
  unit_table <- data.frame(
    unit_id    = c(1L, 2L),
    n_cells    = c(4L, 4L),
    centroid_x = c(1.0, 3.0),
    centroid_y = c(3.0, 1.0)
  )
  grower_type_mat <- matrix(c(1, 1, 0, 0,
                              1, 1, 0, 0,
                              0, 0, 2, 2,
                              0, 0, 2, 2), nrow = nr, byrow = TRUE)
  behavior_params <- list(
    list(detection_prob = 1.0, willingness_to_treat = 1.0,
         treatment_efficacy = 0.8, decision_threshold = 0.0),
    list(detection_prob = 1.0, willingness_to_treat = 1.0,
         treatment_efficacy = 0.8, decision_threshold = 0.0)
  )
  list(
    host_mat         = host_mat,
    infected_mat     = infected_mat,
    unit_mat         = unit_mat,
    unit_table       = unit_table,
    grower_type_mat  = grower_type_mat,
    behavior_params = behavior_params
  )
}

# -- delineate_management_units -----------------------------------------------

test_that("DM1: raster method returns a list with correct elements", {
  r <- make_host_raster()
  res <- delineate_management_units(r, method = "raster")
  expect_type(res, "list")
  expect_named(res, c("management_unit_matrix", "management_unit_table"))
})

test_that("DM2: management_unit_matrix is an integer matrix", {
  r <- make_host_raster()
  res <- delineate_management_units(r, method = "raster")
  expect_true(is.integer(res$management_unit_matrix))
  expect_equal(dim(res$management_unit_matrix), c(4L, 4L))
})

test_that("DM3: raster method produces two distinct units for two disjoint patches", {
  r <- make_host_raster()
  res <- delineate_management_units(r, method = "raster")
  ids <- sort(unique(as.integer(res$management_unit_matrix[res$management_unit_matrix > 0])))
  expect_length(ids, 2L)
})

test_that("DM4: management_unit_table has required columns", {
  r <- make_host_raster()
  res <- delineate_management_units(r, method = "raster")
  expect_true(all(c("unit_id", "n_cells", "centroid_x", "centroid_y") %in%
                    names(res$management_unit_table)))
})

test_that("DM5: grid method produces integer IDs covering host cells", {
  r <- make_host_raster()
  res <- delineate_management_units(r, method = "grid", grid_size = 2)
  expect_true(is.integer(res$management_unit_matrix))
  expect_true(any(res$management_unit_matrix > 0))
})

test_that("DM6: polygon method requires polygon_file", {
  r <- make_host_raster()
  expect_error(
    delineate_management_units(r, method = "polygon", polygon_file = NULL),
    regexp = NULL
  )
})

test_that("DM7: unknown method raises an error", {
  r <- make_host_raster()
  expect_error(delineate_management_units(r, method = "unknown"))
})

test_that("DM8: min_cells removes small patches", {
  r <- make_host_raster()
  res_large <- delineate_management_units(r, method = "raster", min_cells = 1000)
  expect_true(all(res_large$management_unit_matrix == 0L))
})

# -- assign_grower_types -------------------------------------------------------

make_unit_matrix <- function() {
  matrix(c(1, 1, 0, 0,
           1, 1, 0, 0,
           0, 0, 2, 2,
           0, 0, 2, 2), nrow = 4, byrow = TRUE)
}

test_that("GT1: random assignment returns integer matrix of correct dimensions", {
  um <- make_unit_matrix()
  res <- assign_grower_types(
    um, type_names = c("A", "B"), type_probs = c(0.5, 0.5),
    spatial_structure = "random", random_seed = 42
  )
  expect_true(is.integer(res))
  expect_equal(dim(res), dim(um))
})

test_that("GT2: background cells remain 0", {
  um <- make_unit_matrix()
  res <- assign_grower_types(
    um, type_names = c("A", "B"), type_probs = c(0.5, 0.5),
    spatial_structure = "random", random_seed = 1
  )
  expect_true(all(res[um == 0] == 0L))
})

test_that("GT3: assigned types are within 1:K", {
  um <- make_unit_matrix()
  res <- assign_grower_types(
    um, type_names = c("A", "B", "C"), type_probs = c(0.33, 0.33, 0.34),
    spatial_structure = "random", random_seed = 7
  )
  active_vals <- res[um > 0]
  expect_true(all(active_vals %in% 1:3))
})

test_that("GT4: deterministic seed gives reproducible results", {
  um <- make_unit_matrix()
  r1 <- assign_grower_types(um, c("A", "B"), c(0.5, 0.5), random_seed = 99)
  r2 <- assign_grower_types(um, c("A", "B"), c(0.5, 0.5), random_seed = 99)
  expect_identical(r1, r2)
})

test_that("GT5: clustered assignment produces valid type values", {
  um <- make_unit_matrix()
  res <- assign_grower_types(
    um, type_names = c("A", "B"), type_probs = c(0.5, 0.5),
    spatial_structure = "clustered", cluster_range = 2, random_seed = 10
  )
  active <- res[um > 0]
  expect_true(all(active %in% 1:2))
})

test_that("GT6: type_probs must sum to 1", {
  um <- make_unit_matrix()
  expect_error(
    assign_grower_types(um, c("A", "B"), c(0.6, 0.6))
  )
})

# -- grower_decision -----------------------------------------------------------

test_that("GD1: returns a numeric matrix of correct dimensions", {
  s <- make_decision_setup()
  res <- grower_decision(s$infected_mat, s$host_mat, s$grower_type_mat,
                         s$unit_mat, s$unit_table, s$behavior_params,
                         random_seed = 1)
  expect_true(is.matrix(res))
  expect_equal(storage.mode(res), "double")
  expect_equal(dim(res), c(4L, 4L))
})

test_that("GD2: zero infection with high threshold leads to no treatment", {
  s <- make_decision_setup()
  s$infected_mat[] <- 0L
  high_thresh_params <- lapply(s$behavior_params, function(p) {
    p$decision_threshold <- 0.5; p
  })
  res <- grower_decision(s$infected_mat, s$host_mat, s$grower_type_mat,
                         s$unit_mat, s$unit_table, high_thresh_params,
                         random_seed = 1)
  expect_true(all(res == 0))
})

test_that("GD3: binary output contains only 0 and 1", {
  s <- make_decision_setup()
  res <- grower_decision(s$infected_mat, s$host_mat, s$grower_type_mat,
                         s$unit_mat, s$unit_table, s$behavior_params,
                         output = "binary", random_seed = 5)
  expect_true(all(res %in% c(0, 1)))
})

test_that("GD4: continuous output contains only 0 or efficacy values", {
  s <- make_decision_setup()
  res <- grower_decision(s$infected_mat, s$host_mat, s$grower_type_mat,
                         s$unit_mat, s$unit_table, s$behavior_params,
                         output = "continuous", random_seed = 5)
  unique_vals <- unique(as.vector(res))
  expect_true(all(unique_vals %in% c(0, 0.8)))
})

test_that("GD5: treatment only applied to host cells", {
  s <- make_decision_setup()
  res <- grower_decision(s$infected_mat, s$host_mat, s$grower_type_mat,
                         s$unit_mat, s$unit_table, s$behavior_params,
                         random_seed = 1)
  non_host_treated <- res[s$host_mat == 0]
  expect_true(all(non_host_treated == 0))
})

test_that("GD6: invalid output argument raises an error", {
  s <- make_decision_setup()
  expect_error(
    grower_decision(s$infected_mat, s$host_mat, s$grower_type_mat,
                    s$unit_mat, s$unit_table, s$behavior_params,
                    output = "foobar")
  )
})

test_that("GD7: management_unit_table attribute contains infection summary columns", {
  s <- make_decision_setup()
  res <- grower_decision(s$infected_mat, s$host_mat, s$grower_type_mat,
                         s$unit_mat, s$unit_table, s$behavior_params,
                         random_seed = 1)
  tbl <- attr(res, "management_unit_table")
  expect_false(is.null(tbl))
  expect_true(all(c("infected_cells", "host_cells", "prevalence") %in% names(tbl)))
})

test_that("GD8: infected_cells totals match sum of infected matrix", {
  s <- make_decision_setup()
  res <- grower_decision(s$infected_mat, s$host_mat, s$grower_type_mat,
                         s$unit_mat, s$unit_table, s$behavior_params,
                         random_seed = 1)
  tbl <- attr(res, "management_unit_table")
  expect_equal(sum(tbl$infected_cells), sum(s$infected_mat))
})

test_that("GD9: prevalence is between 0 and 1", {
  s <- make_decision_setup()
  res <- grower_decision(s$infected_mat, s$host_mat, s$grower_type_mat,
                         s$unit_mat, s$unit_table, s$behavior_params,
                         random_seed = 1)
  tbl <- attr(res, "management_unit_table")
  expect_true(all(tbl$prevalence >= 0 & tbl$prevalence <= 1))
})

test_that("GD10: reproducible with same random_seed", {
  s <- make_decision_setup()
  r1 <- grower_decision(s$infected_mat, s$host_mat, s$grower_type_mat,
                        s$unit_mat, s$unit_table, s$behavior_params,
                        random_seed = 42)
  r2 <- grower_decision(s$infected_mat, s$host_mat, s$grower_type_mat,
                        s$unit_mat, s$unit_table, s$behavior_params,
                        random_seed = 42)
  expect_identical(r1, r2)
})

# -- behavior_configuration --------------------------------------------------

make_bc_list <- function(...) {
  defaults <- list(
    use_behavior_module     = TRUE,
    grower_id_file           = "grower_id.tif",
    grower_params_file       = "",
    behavior_decision_dates = c("2020-05-01"),
    behavior_params         = list(
      list(detection_prob = 0.8, willingness_to_treat = 0.6,
           treatment_efficacy = 1.0, decision_threshold = 0.05)
    )
  )
  modifyList(defaults, list(...))
}

test_that("BC1: valid config returns a list with no failure", {
  bc <- behavior_configuration(make_bc_list(), "2020-01-01", "2020-12-31")
  expect_null(bc$failure)
})

test_that("BC2: use_behavior_module defaults to FALSE when absent", {
  bc_in <- make_bc_list()
  bc_in$use_behavior_module <- NULL
  bc <- behavior_configuration(bc_in, "2020-01-01", "2020-12-31")
  expect_false(bc$use_behavior_module)
})

test_that("BC3: legacy enable_behavior_module is accepted as alias", {
  bc_in <- make_bc_list()
  bc_in$use_behavior_module <- NULL
  bc_in$enable_behavior_module <- TRUE
  bc <- behavior_configuration(bc_in, "2020-01-01", "2020-12-31")
  expect_true(bc$use_behavior_module)
})

test_that("BC4: decision date outside simulation window returns failure", {
  bc_in <- make_bc_list(behavior_decision_dates = c("2025-01-01"))
  bc <- behavior_configuration(bc_in, "2020-01-01", "2020-12-31")
  expect_false(is.null(bc$failure))
})

test_that("BC5: grower_id_file and grower_params_file are preserved", {
  bc_in <- make_bc_list(grower_id_file = "my_growers.tif",
                        grower_params_file = "my_params.csv")
  bc <- behavior_configuration(bc_in, "2020-01-01", "2020-12-31")
  expect_equal(bc$grower_id_file, "my_growers.tif")
  expect_equal(bc$grower_params_file, "my_params.csv")
})

test_that("BC6: missing behavior_params sub-fields receive defaults", {
  bc_in <- make_bc_list(behavior_params = list(list()))
  bc <- behavior_configuration(bc_in, "2020-01-01", "2020-12-31")
  expect_equal(bc$behavior_params[[1]]$detection_prob,       1.0)
  expect_equal(bc$behavior_params[[1]]$willingness_to_treat, 0.5)
  expect_equal(bc$behavior_params[[1]]$treatment_efficacy,   1.0)
  expect_equal(bc$behavior_params[[1]]$decision_threshold,   0.0)
})

test_that("BC7: legacy decision_dates field is accepted as alias", {
  bc_in <- make_bc_list()
  bc_in$behavior_decision_dates <- NULL
  bc_in$decision_dates <- c("2020-06-01")
  bc <- behavior_configuration(bc_in, "2020-01-01", "2020-12-31")
  expect_equal(bc$behavior_decision_dates, c("2020-06-01"))
})
