# These functions are designed to reduce code complexity and the need to copy
# and past code across main functions

# calibration success metric option
success_metric_options <- c("quantity", "allocation", "configuration", "quantity and allocation",
                            "quantity and configuration", "allocation and configuration",
                            "quantity, allocation, and configuration", "accuracy", "precision",
                            "recall", "specificity", "accuracy and precision",
                            "accuracy and specificity", "accuracy and recall",
                            "precision and recall", "precision and specificity",
                            "recall and specificity", "accuracy, precision, and recall",
                            "accuracy, precision, and specificity",
                            "accuracy, recall, and specificity",
                            "precision, recall, and specificity",
                            "accuracy, precision, recall, and specificity",
                            "rmse", "distance", "mcc", "mcc and quantity", "mcc and distance",
                            "rmse and distance", "mcc and configuration", "mcc and RMSE",
                            "mcc, quantity, and configuration")

quantity_list <- c("quantity", "quantity and allocation", "quantity and configuration",
                   "quantity, allocation, and configuration", "mcc and quantity",
                   "mcc, quantity, and configuration")

allocation_list <- c("allocation", "quantity and allocation", "allocation and configuration",
                     "quantity, allocation, and configuration")

configuration_list <- c("configuration",  "quantity and configuration",
                        "allocation and configuration", "quantity, allocation, and configuration",
                        "mcc and configuration", "mcc, quantity, and configuration")

accurracy_list <- c("accuracy", "accuracy and precision", "accuracy and specificity",
                    "accuracy and recall", "accuracy, precision, and recall",
                    "accuracy, precision, and specificity",
                    "accuracy, recall, and specificity",
                    "accuracy, precision, recall, and specificity")

precision_list <- c("precision", "accuracy and precision", "precision and recall",
                    "precision and specificity", "accuracy, precision, and recall",
                    "accuracy, precision, and specificity",
                    "precision, recall, and specificity",
                    "accuracy, precision, recall, and specificity")

recall_list <- c("recall", "accuracy and recall", "precision and recall", "recall and specificity",
                 "accuracy, precision, and recall", "accuracy, recall, and specificity",
                 "precision, recall, and specificity",
                 "accuracy, precision, recall, and specificity")

specificity_list <- c("specificity", "accuracy and specificity", "precision and specificity",
                      "recall and specificity", "accuracy, precision, and specificity",
                      "accuracy, recall, and specificity",  "precision, recall, and specificity",
                      "accuracy, precision, recall, and specificity")

rmse_list <- c("rmse", "rmse and distance", "mcc and RMSE")

distance_list <- c("distance", "mcc and distance", "rmse and distance")

mcc_list <- c("mcc", "mcc and quantity", "mcc and distance", "mcc and configuration",
              "mcc and RMSE", "mcc, quantity, and configuration")

weather_type_list <- c("deterministic", "probabilistic", "none")

"%notin%" <- Negate("%in%")

set_success_metrics <- function(config) {
  config$use_quantity <- FALSE
  config$use_allocation <- FALSE
  config$use_configuration <- FALSE
  config$use_accuracy <- FALSE
  config$use_precision <- FALSE
  config$use_recall <- FALSE
  config$use_specificity <- FALSE
  config$use_rmse <- FALSE
  config$use_distance <- FALSE
  config$use_mcc <- FALSE

  if (config$success_metric %in% quantity_list) {
    config$use_quantity <- TRUE
  }

  if (config$success_metric %in% allocation_list) {
    config$use_allocation <- TRUE
  }

  if (config$success_metric %in% configuration_list) {
    config$use_configuration <- TRUE
  }

  if (config$success_metric %in% accurracy_list) {
    config$use_accuracy <- TRUE
  }

  if (config$success_metric %in% precision_list) {
    config$use_precision <- TRUE
  }

  if (config$success_metric %in% recall_list) {
    config$use_recall <- TRUE
  }

  if (config$success_metric %in% specificity_list) {
    config$use_specificity <- TRUE
  }

  if (config$success_metric %in% rmse_list) {
    config$use_rmse <- TRUE
  }

  if (config$success_metric %in% distance_list) {
    config$use_distance <- TRUE
  }

  if (config$success_metric %in% mcc_list) {
    config$use_mcc <- TRUE
  }

  return(config)
}

# creates a matrix from a matrix of mean values and a matrix of standard deviations. The two
# matrices must be the same size.
matrix_norm_distribution <- function(mean_matrix, sd_matrix) {
  new_matrix <-
    round(matrix(mapply(function(x, y) {rnorm(x, y, n = 1)}, x = mean_matrix, y = sd_matrix),
                 nrow(mean_matrix), ncol(mean_matrix)))
  new_matrix[is.na(new_matrix)] <- 0
  new_matrix[new_matrix < 0] <- 0
  return(new_matrix)
}

# Uncertainty propagation for raster data sets, expects a spatRaster with 2
# layers (mean and standard deviation)
output_from_raster_mean_and_sd <- function(x) {
  x[[1]] <- terra::classify(x[[1]], matrix(c(-Inf, 0, 0), ncol = 3, byrow = TRUE))
  x[[2]] <- terra::classify(x[[2]], matrix(c(-Inf, 0, 0), ncol = 3, byrow = TRUE))
  fun <- function(x) {
    round(rnorm(1, mean = x[1], sd = x[2]), digits = 0)
  }
  x2 <- suppressWarnings(terra::app(x, fun))
  return(x2)
}

# function for getting all infected locations based on rook or queens rule for
# assessing clusters of infections.
get_all_infected <- function(rast, direction = 4) {
  # get infections as points
  p <- terra::as.points(rast)
  rast <- terra::classify(rast, matrix(c(0, NA), ncol = 2, byrow = TRUE), right = NA)
  names(rast) <- "group"
  names(p) <- "data"
  p <- p[p$data > 0]
  infections <- data.frame(terra::extract(rast, p, cells = TRUE))
  infections <- infections[, 2:3]
  names(infections) <- c("detections", "cells")
  if (direction %in% c(4, 8)) {
    infections$i <- terra::colFromCell(rast, infections$cells)
    infections$j <- terra::rowFromCell(rast, infections$cells)
    r <- terra::patches(rast, direction = direction, zeroAsNA = TRUE)
    infections$group <- terra::extract(r, p)$patches
  } else {
    return("direction should be either of 4 or 8")
  }

  infections$group_size <- 0
  groups <- data.frame(table(infections$group))
  groups$Var1 <- as.numeric(groups$Var1)
  for (m in seq_len(nrow(groups))) {
    infections$group_size[infections$group == groups$Var1[m]] <- groups$Freq[m]
  }
  names(infections) <-
    c("detections", "cell_number", "i", "j", "group", "group_size")
  return(infections)
}

# returns the foci of infestation for a spatRaster Object
get_foci <- function(rast) {
  indexes <- get_all_infected(rast)
  center <-
    data.frame(i = round(mean(indexes$i), digits = 0), j = round(mean(indexes$j, digits = 0)))
  return(center)
}

# returns the border of the infected area for a spatRaster
get_infection_border <- function(rast) {
  s <- get_all_infected(rast)
  min_max_col <-
    data.frame(row_number = seq(1, nrow(rast), 1),
               min_j = rep(0, nrow(rast)),
               max_j = rep(0, nrow(rast)))
  min_max_row <-
    data.frame(column_number = seq(1, ncol(rast), 1),
               min_i = rep(0, ncol(rast)),
               max_i = rep(0, ncol(rast)))
  for (i in seq_len(nrow(rast))) {
    if (length(s$j[s$i %in% i]) > 0) {
      min_max_col$max_j[i] <- max(s$j[s$i %in% i])
      min_max_col$min_j[i] <- min(s$j[s$i %in% i])
    } else {
      min_max_col$max_j[i] <- 0
      min_max_col$min_j[i] <- 0
    }
  }

  for (j in seq_len(ncol(rast))) {
    if (length(s$i[s$j %in% j]) > 0) {
      min_max_row$max_i[j] <- max(s$i[s$j %in% j])
      min_max_row$min_i[j] <- min(s$i[s$j %in% j])
    } else {
      min_max_row$max_i[j] <- 0
      min_max_row$min_i[j] <- 0
    }
  }

  border1 <-
    data.frame(rows = rep(min_max_col$row_number, 2),
               cols = c(min_max_col$min_j, min_max_col$max_j))
  border2 <-
    data.frame(rows = c(min_max_row$min_i, min_max_row$max_i),
               cols = rep(min_max_row$column_number))
  border <- rbind(border1, border2)
  borders <- data.frame(table(border[, 1:2]))
  borders <- borders[borders$Freq > 0, ]
  borders <- borders[borders$rows != 0, ]
  borders <- borders[borders$cols != 0, ]
  borders$rows <- as.numeric(as.character(borders$rows))
  borders$cols <- as.numeric(as.character(borders$cols))

  names(borders) <- c("i", "j", "freq")
  return(borders)
}

# returns the distances of points from either the Foci, the Border, or a set of
# points
get_infection_distances <- function(rast, method = "Foci", points = c()) {
  infections <- get_all_infected(rast)
  if (method == "Foci") {
    points <- get_foci(rast)
  } else if (method == "Border") {
    points <- get_infection_border(rast)
  } else if (method == "Points") {
    points <- points
  } else {
    return("method must be one of 'Foci', 'Border', or'Points'")
  }
  for (i in seq_len(nrow(infections))) {
    minimum_distance <- Inf
    for (l in seq_len(nrow(points))) {
      new_distance <-
        sqrt((abs(infections$i[i] - points$i[l]))^2 + (abs(infections$j[i] - points$j[l]))^2)
      if (new_distance < minimum_distance) {
        minimum_distance <- new_distance
      }
    }
    infections$distance[i] <- minimum_distance
  }

  return(infections)
}

# returns a set of treatments for a group of species infestations based on
# multiple rules
treatment_auto <- function(rasts,
                           rasts2,
                           method = "Foci",
                           priority = "group size",
                           number_of_locations = 1,
                           points = c(),
                           treatment_efficacy = 1,
                           buffer_cells = 1.5,
                           direction_first = TRUE,
                           treatment_priority = "equal",
                           treatment_rank = c(0)) {
  # get distances and groups and group size
  if (method == "Points") {
    points <- points
  }

  if (treatment_priority == "equal") {
    rasts <-
      terra::tapp(rasts,
                  index = rep(1, terra::nlyr(rasts)),
                  fun = sum)
    rasts2 <-
      terra::tapp(rasts2,
                  index = rep(1, terra::nlyr(rasts2)),
                  fun = sum)
  } else if (treatment_priority == "ranked") {
    for (r in seq_len(length(treatment_rank))) {
      if (r == 1) {
        raste <- rast(rasts[[match(r, treatment_rank)]])
        terra::values(raste) <- terra::values(rasts[[match(r, treatment_rank)]])
        raste2 <- rast(rasts2[[match(r, treatment_rank)]])
        terra::values(raste2) <-
          terra::values(rasts2[[match(r, treatment_rank)]])
      } else if (r > 1) {
        raste <- c(raste, rasts[[match(r, treatment_rank)]])
        raste2 <- c(raste2, rasts2[[match(r, treatment_rank)]])
      }
    }
    rasts <- raste
    rasts2 <- raste2
  }

  total_infs <- c(0)
  cells_treated <- 0
  treatment <- rasts[[1]]
  treatment[] <- 0
  names(treatment) <- treatment
  for (q in 1:terra::nlyr(rasts)) {
    rast <- rasts[[q]]
    rast2 <- rasts2[[q]]
    total_infs[q] <- sum(rast[rast > 0] > 0)
    if (q > 1) {
      treatment <- treatment
      cells_treated <- cells_treated
    }

    if (total_infs[q] > 0 && cells_treated[[1]] < number_of_locations) {

      infections <- get_infection_distances(rast = rast, method = method, points = points)
      infections$host <- 0
      for (p in seq_len(nrow(infections))) {
        infections$host[p] <- rast2[infections$i[p], infections$j[p]]
      }

      if (priority == "group size") {
        if (direction_first) {
          infections <-
            infections[order(infections$distance,
                             infections$group_size,
                             decreasing = c(FALSE, TRUE)), ]
        } else {
          infections <-
            infections[order(infections$group_size,
                             infections$distance,
                             decreasing = c(TRUE, FALSE)), ]
        }
        groups_unique <- unique(infections$group)
        for (group in groups_unique) {
          managed_group <- infections[infections$group == group, ]
          for (m in seq_len(nrow(managed_group))) {
            i <- managed_group$i[m]
            j <- managed_group$j[m]
            if (treatment[i, j] < 1 & (rast[i, j] | rast2[i, j])) {
              value <- min(1, treatment[i, j]$treatment + 1)
              if (value > treatment[i, j]) {
                cells_treated <- cells_treated + value - treatment[i, j]
              } else {
                cells_treated <- cells_treated + value
              }
              treatment[i, j] <- value
              if (cells_treated >= number_of_locations) {
                break
                }
            }
          }
          for (m in seq_len(nrow(managed_group))) {
            i <- managed_group$i[m]
            j <- managed_group$j[m]
            i_s <- seq(floor(i - buffer_cells), ceiling(i + buffer_cells), 1)
            j_s <- seq(floor(j - buffer_cells), ceiling(j + buffer_cells), 1)
            i_s <- i_s[i_s > 0]
            j_s <- j_s[j_s > 0]
            for (s in seq_len(length(i_s))) {
              for (n in seq_len(length(j_s))) {
                if (treatment[i_s[s], j_s[n]] < 1 &
                    (rast[i_s[s], j_s[n]] | rast2[i_s[s], j_s[n]])) {
                  if (cells_treated >= number_of_locations) {
                    break
                  }
                  if (abs(i - i_s[s]) > buffer_cells |
                      abs(j - j_s[n]) > buffer_cells) {
                    value <-
                      min(1, treatment[i_s[s], j_s[n]]$treatment +
                            (buffer_cells - floor(buffer_cells)))
                    if (value > treatment[i_s[s], j_s[n]]) {
                      cells_treated <- cells_treated + value - treatment[i_s[s], j_s[n]]
                    } else {
                      cells_treated <- cells_treated + value
                    }
                    treatment[i_s[s], j_s[n]] <- value
                  } else if (rast[i_s[s], j_s[n]] | rast2[i_s[s], j_s[n]]) {
                    if (treatment[i_s[s], j_s[n]] < 1) {
                      if (cells_treated >= number_of_locations) {
                        break
                      }
                      cells_treated <- cells_treated + (1 - treatment[i_s[s], j_s[n]])
                      treatment[i_s[s], j_s[n]] <- 1
                    }
                  }
                }
                if (cells_treated >= number_of_locations) {
                  break
                  }
              }
              if (cells_treated >= number_of_locations) {
                break
                }
            }
            if (cells_treated >= number_of_locations) {
              break
              }
          }
          if (cells_treated >= number_of_locations) {
            break
            }
        }
      } else if (priority == "host") {
        if (direction_first) {
          infections <-
            infections[order(infections$distance,
                             infections$host,
                             decreasing = c(FALSE, TRUE)), ]
        } else {
          infections <-
            infections[order(infections$host,
                             infections$distance,
                             decreasing = c(TRUE, FALSE)), ]
        }
        for (t in seq_len(nrow(infections))) {
          i <- infections$i[t]
          j <- infections$j[t]
          if (treatment[i, j] < 1) {
            cells_treated <- cells_treated + (1 - treatment[i, j])
            treatment[i, j] <- 1
          }
          i_s <- seq(floor(i - buffer_cells), ceiling(i + buffer_cells), 1)
          j_s <- seq(floor(j - buffer_cells), ceiling(j + buffer_cells), 1)
          i_s <- i_s[i_s > 0]
          j_s <- j_s[j_s > 0]
          for (s in seq_len(length(i_s))) {
            for (n in seq_len(length(j_s))) {
              if (cells_treated >= number_of_locations) {
                break
                }
              if (treatment[i_s[s], j_s[n]] < 1 &
                  (rast[i_s[s], j_s[n]] | rast2[i_s[s], j_s[n]])) {
                if (abs(i - i_s[s]) > buffer_cells |
                    abs(j - j_s[n]) > buffer_cells) {
                  value <-
                    min(1, treatment[i_s[s], j_s[n]]$treatment +
                          (buffer_cells - floor(buffer_cells)))
                  if (value > treatment[i_s[s], j_s[n]]) {
                    cells_treated <- cells_treated + value - treatment[i_s[s], j_s[n]]
                  } else {
                    cells_treated <- cells_treated + value
                  }
                  treatment[i_s[s], j_s[n]] <- value
                } else if (rast[i_s[s], j_s[n]] | rast2[i_s[s], j_s[n]]) {
                  if (treatment[i_s[s], j_s[n]] < 1) {
                    cells_treated <- cells_treated + (1 - treatment[i_s[s], j_s[n]])
                    treatment[i_s[s], j_s[n]] <- 1
                  }
                }
              }
            }
          }
          if (cells_treated >= number_of_locations) {
            break
          }
        }
      } else if (priority == "infected") {
        if (direction_first) {
          infections <-
            infections[order(infections$distance,
                             infections$detections,
                             decreasing = c(FALSE, TRUE)), ]
        } else {
          infections <-
            infections[order(infections$detections,
                             infections$distance,
                             decreasing = c(TRUE, FALSE)), ]
        }
        for (t in seq_len(nrow(infections))) {
          i <- infections$i[t]
          j <- infections$j[t]
          if (treatment[i, j] < 1) {
            cells_treated <- cells_treated + (1 - treatment[i, j])
            treatment[i, j] <- 1
          }
          i_s <- seq(floor(i - buffer_cells), ceiling(i + buffer_cells), 1)
          j_s <- seq(floor(j - buffer_cells), ceiling(j + buffer_cells), 1)
          i_s <- i_s[i_s > 0]
          j_s <- j_s[j_s > 0]
          for (s in seq_len(length(i_s))) {
            for (n in seq_len(length(j_s))) {
              if (cells_treated >= number_of_locations) {
                break
              }
              if (treatment[i_s[s], j_s[n]] < 1 &
                  (rast[i_s[s], j_s[n]] | rast2[i_s[s], j_s[n]])) {
                if (abs(i - i_s[s]) > buffer_cells |
                    abs(j - j_s[n]) > buffer_cells) {
                  value <-
                    min(1, treatment[i_s[s], j_s[n]]$treatment +
                          (buffer_cells - floor(buffer_cells)))
                  if (value > treatment[i_s[s], j_s[n]]) {
                    cells_treated <- cells_treated + value - treatment[i_s[s], j_s[n]]
                  } else {
                    cells_treated <- cells_treated + value
                  }
                  treatment[i_s[s], j_s[n]] <- value
                } else if (rast[i_s[s], j_s[n]] | rast2[i_s[s], j_s[n]]) {
                  if (treatment[i_s[s], j_s[n]] < 1) {
                    cells_treated <- cells_treated + (1 - treatment[i_s[s], j_s[n]])
                    treatment[i_s[s], j_s[n]] <- 1
                  }
                }
              }

            }
          }
          if (cells_treated >= number_of_locations) {
            break
          }
        }
      } else {
        return('priority needs to be one of "group size", "host", or "infected"')
      }
    }
  }
  treatment <- treatment[[1]]

  treatment <- treatment * treatment_efficacy

  return(treatment)
}
