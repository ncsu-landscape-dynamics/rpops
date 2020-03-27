## These functions are designed to reduce code complexity and the need to copy and past code across main functions

getAllInfected <- function(rast, direction = 4) {
  p <- rasterToPoints(rast, fun=function(x){x>0}, spatial = TRUE)
  infections <- data.frame(extract(rast, p, cellnumbers=TRUE))
  if (direction == 4) {
    cellnumbersb <- data.frame(extract(rast, p, buf = xres(rast), cellnumbers=TRUE))
  } else if (direction == 8) {
    cellnumbersb <- data.frame(extract(rast, p, buf = 1.5 * xres(rast), cellnumbers=TRUE))
  } else {
    return("direction should be either of 4 or 8")
  }
  infections$i <- infections$cells %/% ncol(rast) + 1
  infections$j <- infections$cells %% ncol(rast)
  infections$group <- seq(1, nrow(infections),1)
  cell_cols <- seq(1, ncol(cellnumbersb)-1,2)
  for (i in cell_cols) {
    infections$group[infections$cells %in% cellnumbersb[,i]] <- min(infections$group[infections$cells %in% cellnumbersb[,i]])
    
  }
  groups <- unique(infections$group)
  group_map <- seq(1, length(groups), 1)
  for (m in 1:length(groups)) {
    infections$group_size[infections$group == groups[m]] <- nrow(infections[infections$group == groups[m],])
    infections$group[infections$group == groups[m]] <- group_map[m]
    
  }
  names(infections) <- c('cell_number', 'detections', 'i', 'j', 'group', 'group_size')
  return(infections)
}

getFoci <- function(rast) {
  indexes <- getAllInfected(rast)
  center <- data.frame(i = floor(mean(indexes$i)), j = floor(mean(indexes$j)))
  return(center)
}

getInfectionBorder <- function(rast) {
  s <- getAllInfected(rast)
  min_max_col <- data.frame(row_number = seq(1, nrow(rast),1), min_j = rep(0,nrow(rast)), max_j = rep(0,nrow(rast)))
  min_max_row <- data.frame(column_number = seq(1, ncol(rast), 1), min_i = rep(0,ncol(rast)), max_i = rep(0,ncol(rast)))
  for (i in 1:nrow(rast)) {
    if (length(s$j[s$i %in% i]) > 0) {
      min_max_col$max_j[i] <- max(s$j[s$i %in% i])
      min_max_col$min_j[i] <- min(s$j[s$i %in% i])
    } else {
      min_max_col$max_j[i] <- 0
      min_max_col$min_j[i] <- 0
    }
  }
  
  for (j in 1:ncol(rast)) {
    if (length(s$i[s$j %in% j]) > 0) {
      min_max_row$max_i[j] <- max(s$i[s$j %in% j])
      min_max_row$min_i[j] <- min(s$i[s$j %in% j])
    } else {
      min_max_row$max_i[j] <- 0
      min_max_row$min_i[j] <- 0
    }
  }
  
  border1 <- data.frame(rows = rep(min_max_col$row_number,2), cols = c(min_max_col$min_j, min_max_col$max_j))
  border2 <- data.frame(rows = c(min_max_row$min_i, min_max_row$max_i), cols = rep(min_max_row$column_number))
  border <- rbind(border1, border2)
  borders <- data.frame(table(border[,1:2]))
  borders <- borders[borders$Freq > 0, ]
  borders <- borders[borders$rows != 0, ]
  borders <- borders[borders$cols != 0, ]
  borders$rows <- as.numeric(as.character(borders$rows))
  borders$cols <- as.numeric(as.character(borders$cols))
  
  names(borders) <- c('i','j', 'freq')
  return(borders)
}

## options for get 
getInfectionDistances <- function(rast, method = 'Foci', points = c()) {
  infections <- getAllInfected(rast)
  if (method == "Foci") {
    points <- getFoci(rast)
  } else if (method == "Border") {
    points <- getInfectionBorder(rast)
  } else if (method == "Points") {
    points <- points
  } else {
    return("method must be one of 'Foci', 'Border', or'Points'")
  }
  # minimum_distance <- Inf
  for (i in 1:nrow(infections)) {
    minimum_distance <- Inf
    for (l in 1:nrow(points)) {
      new_distance <- sqrt((abs(infections$i[i] - points$i[l]))^2 + (abs(infections$j[i] - points$j[l]))^2)
      if (new_distance < minimum_distance) {
        minimum_distance <- new_distance
      }
    }
    infections$distance[i] <- minimum_distance
  }
  
  return(infections)
}

treatmentAuto <- function(rasts, rasts2, method = 'Foci', priority = 'group size', 
                          number_of_locations = 1, points = c(), treatment_efficacy = 1, 
                          buffer_cells = 1.5, direction_first = TRUE, treatment_priority = "equal",
                          treatment_rank = c(0)) {
  ## get distances and groups and group size
  if (method == "Points") {
    points <- points
  }
  
  if (treatment_priority == "equal") {
    rasts <- raster::stackApply(rasts, indices = rep(1, raster::nlayers(rasts)), fun = sum)
    rasts2 <- raster::stackApply(rasts2, indices = rep(1, raster::nlayers(rasts2)), fun = sum)
  } else if (treatment_priority == "ranked") {
    if (all(treatment_rank == c(1, 0))) {
      rasts <- rasts
      rasts2 <- rasts2
    } else if (all(treatment_rank == c(0, 1))) {
      rasts <- stack(rasts[[2]], rasts[[1]])
      rasts2 <- stack(rasts2[[2]], rasts2[[1]])
    }
  }

  total_infs <- c(0)
  cells_treated <- 0
  treatments <- data.frame(i = 0, j = 0, value = 0)
  treatment <- rasts
  treatment[] <- 0
  for (q in 1:nlayers(rasts)) {
    rast <- rasts[[q]]
    rast2 <- rasts2[[q]]
    total_infs[q] <- sum(rast[rast > 0] > 0)
    if (q > 1) {
      treatment <- treatment
      cells_treated <- cells_treated
    }
    
    if (total_infs[q] > 0 && cells_treated < number_of_locations){
      
      infections <- getInfectionDistances(rast = rast, method = method, points = points)
      infections$host <- 0
      for (p in 1:length(infections)) {
        infections$host <- rast2[infections$i[p], infections$j[p]]
      }
      
      if (priority == 'group size') {
        if (direction_first) {
          infections <- infections[order(infections$distance, infections$group_size, decreasing = c(FALSE, TRUE)),]
        } else {
          infections <- infections[order(infections$group_size, infections$distance, decreasing = c(TRUE, FALSE)),]
        }
        groups_unique <- unique(infections$group)
        for (group in groups_unique) {
          managed_group <- infections[infections$group == group,]
          group_size <- nrow(managed_group)
          for (m in 1:nrow(managed_group)) {
            i = managed_group$i[m]
            j = managed_group$j[m]
            if (treatment[i,j] < 1 & (rast[i, j] | rast2[i, j])) {
              value <- min(1, treatment[i, j] + 1)
              if (value > treatment[i, j]) {
                cells_treated <- cells_treated + value - treatment[i, j]
              } else {
                cells_treated <- cells_treated + value
              }
              treatment[i, j] <- value
              if (cells_treated >= number_of_locations) {break}
            } 
          }
          for (m in 1:nrow(managed_group)) {
            i = managed_group$i[m]
            j = managed_group$j[m]
            i_s <- seq(floor(i - buffer_cells), ceiling(i + buffer_cells), 1)
            j_s <- seq(floor(j - buffer_cells), ceiling(j + buffer_cells), 1)
            for (s in 1:length(i_s)) {
              for(n in 1:length(j_s)) {
                if (treatment[i_s[s], j_s[n]] < 1 & (rast[i_s[s], j_s[n]] | rast2[i_s[s], j_s[n]])) {
                  if (abs(i - i_s[s]) > buffer_cells | abs(j - j_s[n]) > buffer_cells) {
                    value <- min(1, treatment[i_s[s], j_s[n]] + (buffer_cells - floor(buffer_cells)))
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
                  if (cells_treated >= number_of_locations) {break}
                }
                if (cells_treated >= number_of_locations) {break}
              }
              if (cells_treated >= number_of_locations) {break}
            }
            if (cells_treated >= number_of_locations) {break}
          }
          if (cells_treated >= number_of_locations) {break}
        }
      } else if (priority == 'host') {
        if (direction_first) {
          infections <- infections[order(infections$distance, infections$host, decreasing = c(FALSE, TRUE)),]
        } else {
          infections <- infections[order(infections$host, infections$distance, decreasing = c(TRUE, FALSE)),]
        }
        for (t in 1:nrow(infections)) {
          i <- infections$i[t]
          j <- infections$j[t]
          if (treatment[i, j] < 1) {
            cells_treated <- cells_treated + (1 - treatment[i, j])
            treatment[i, j] <- 1
          }
          i_s <- seq(floor(i - buffer_cells), ceiling(i + buffer_cells), 1)
          j_s <- seq(floor(j - buffer_cells), ceiling(j + buffer_cells), 1)
          for (s in 1:length(i_s)) {
            for(n in 1:length(j_s)) {
              if (cells_treated >= number_of_locations) {break}
              if (treatment[i_s[s], j_s[n]] < 1 & (rast[i_s[s], j_s[n]] | rast2[i_s[s], j_s[n]])) {
                if (abs(i - i_s[s]) > buffer_cells | abs(j - j_s[n]) > buffer_cells) {
                  value <- min(1, treatment[i_s[s], j_s[n]] + (buffer_cells - floor(buffer_cells)))
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
          if (cells_treated >= number_of_locations) {break}
        }
      } else if (priority == 'infected') {
        if (direction_first) {
          infections <- infections[order(infections$distance, infections$detections, decreasing = c(FALSE, TRUE)),]
        } else {
          infections <- infections[order(infections$detections, infections$distance, decreasing = c(TRUE, FALSE)),]
        }
        for (t in 1:nrow(infections)) {
          i <- infections$i[t]
          j <- infections$j[t]
          if (treatment[i, j] < 1) {
            cells_treated <- cells_treated + (1 - treatment[i, j])
            treatment[i, j] <- 1
          }
          i_s <- seq(floor(i - buffer_cells), ceiling(i + buffer_cells), 1)
          j_s <- seq(floor(j - buffer_cells), ceiling(j + buffer_cells), 1)
          for (s in 1:length(i_s)) {
            for(n in 1:length(j_s)) {
              if (cells_treated >= number_of_locations) {break}
              if (treatment[i_s[s], j_s[n]] < 1 & (rast[i_s[s], j_s[n]] | rast2[i_s[s], j_s[n]])) {
                if (abs(i - i_s[s]) > buffer_cells | abs(j - j_s[n]) > buffer_cells) {
                  value <- min(1, treatment[i_s[s], j_s[n]] + (buffer_cells - floor(buffer_cells)))
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
          if (cells_treated >= number_of_locations) {break}
        }
      } else {
        return('priority needs to be one of "group size", "host", or "infected"')
      }
    }
    print(q)
  }
  treatment <- treatment[[1]]

  treatment <- treatment * treatment_efficacy
  
  return(treatment)
}
