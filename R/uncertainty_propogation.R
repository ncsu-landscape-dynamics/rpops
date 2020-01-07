## Uncertainty propogation checks and passing

output_from_raster_mean_and_sd <- function(x) {
  x[[1]] <- reclassify(x[[1]], matrix(c(-Inf, 0, 0), ncol = 3, byrow = TRUE))
  x[[2]] <- reclassify(x[[2]], matrix(c(-Inf, 0, 0), ncol = 3, byrow = TRUE))
  # x <- reclassify(x, c(-Inf, 0, 0))
  fun <- function(x) {round(rnorm(1, mean = x[1], sd = x[2]), digits = 0)}
  x2 <- suppressWarnings(calc(x, fun))
  return(x2)
}

choose_variable_based_on_probability <- function(x, n = 1) {
  names(x) <- c('var', 'prob')
  for (i in 1:nrow(x)){
    x$bin[i] <- sum(x$prob[1:i])
    if (i == nrow(x)) {
      if ((x$bin[i] < 1.00 && x$bin[i] > 0.985) || (x$bin[i] > 1.00 && x$bin[i] < 1.01)) {
        x$bin[i] <- 1.00
      }
    }
  }
    rn <- runif(n)
    value <- rep(0,n)
    for (j in 1:length(rn)){
      for (i in 1:nrow(x)) {
        if (i <= 1) {
          if (rn[j] <= x$bin[i]) {
            value[j] <- x$var[i]
          }
        } else {
          if (rn[j] <= x$bin[i] && rn > x$bin[i-1]) {
            value[j] <- x$var[i]
          }
        }
      }
    }
  return(value)
}

uncertainty_check <- function(priors, round_to = 0, n = 1) {
  checks_passed <- TRUE
  
  if (class(priors) == "numeric" && length(priors) == 2) {
    priors <- matrix(priors, ncol = 2)
  } 
  
  if (class(priors) %in% c("matrix", "data.frame") && ncol(priors) == 2) {
    if (class(priors) == "matrix" && nrow(priors) == 1) {
      value <- round(rnorm(n, priors[1], priors[2]), digits = round_to)
    } else if (class(priors) %in% c("matrix", "data.frame") && nrow(priors) >= 1) {
      if (class(priors) == "matrix") {
        priors <- data.frame(priors)
      }
      if (round(sum(priors[,2]), 3) < 0.985 || round(sum(priors[,2]), 3) > 1.01) {
        checks_passed <- FALSE
        failed_check <- "probabilities do not add up to 1"
      }
      value <- choose_variable_based_on_probability(priors, n)
    }
  } else if (class(priors) == "numeric" && length(priors) == 1) {
    value <- rep(priors, n)
  } else {
    checks_passed <- FALSE
    failed_check <- "Incorrect format for priors"
  }
  
  if (checks_passed) {
    outs <- list(checks_passed, value)
    names(outs) <- c('checks_passed', 'value')
    return(outs)
  } else {
    outs <- list(checks_passed, failed_check)
    names(outs) <- c('checks_passed', 'failed_check')
    return(outs)
  }
}

