###############################################################################
# Weighted and un-weighted MA and RMA functions and corresponding plots for Gačnik et al. (2025)
# CITATION: Gačnik J., Žagar K., Hatvani I. G., Kern Z., and Vreča P. (in preparation): Climate change reflected in 40-year isotopic composition trends of precipitation in Slovenia
# CODE AUTHOR: JAN GAČNIK, MAY 2025
# R version: 4.4.2
###############################################################################
# Libraries used for the functions below are:
invisible(lapply(c("ggplot2", "dplyr"), library, character.only = TRUE))

# Calculate MA and weighted MA regressions
# For the RMA_weighted function, x, y, and weights need to be filtered for NAs
RMA_weighted <- function(x, y, weights = NULL) {
  if (!is.null(weights)) {
    weight_factors <- weights/sum(weights)
    n <- length(x)
    x_mean <- mean(x)
    y_mean <- mean(y)
    x_mean_weighted <- sum(x * weight_factors)
    y_mean_weighted <- sum(y * weight_factors)
    sum_x <- sum(x)
    sum_y <- sum(y)
    sum_x2 <- sum(x^2)
    sum_y2 <- sum(y^2)
    sum_u2 <- sum((x - x_mean)^2)
    sum_v2 <- sum((y - y_mean)^2)
    sum_uv <- sum((x - x_mean) * (y - y_mean))
    u_w <- x - x_mean_weighted
    v_w <- y - y_mean_weighted
    slope_RMA_normal <- sqrt((sum_y2 - sum_y^2/n) / (sum_x2 - sum_x^2/n))
    intercept_RMA_normal <- y_mean - slope_RMA_normal * x_mean
    error_slope_RMA_normal <- sqrt( ((sum((y - (slope_RMA_normal * x + intercept_RMA_normal))^2))/(n - 2)) / sum((x - x_mean)^2) )
    error_intercept_RMA_normal <- error_slope_RMA_normal * sqrt(sum_x2 / n)
    slope_MA_normal <- (sum_v2 - sum_u2 + sqrt((sum_v2 - sum_u2)^2 + 4 * ((sum_uv)^2))) / (2 * sum_uv)
    intercept_MA_normal <- y_mean - slope_MA_normal * x_mean
    error_slope_MA_normal <- sqrt( ((sum((y - (slope_MA_normal * x + intercept_MA_normal))^2))/(n - 2)) / sum((x - x_mean)^2) )
    error_intercept_MA_normal <- error_slope_MA_normal * sqrt(sum_x2 / n)
    slope_RMA_weighted <- sqrt(sum(weight_factors * v_w^2) / sum(weight_factors * u_w^2))
    intercept_RMA_weighted <- y_mean_weighted - slope_RMA_weighted * x_mean_weighted
    error_slope_RMA_weighted <- sqrt( ((n * sum(weights * ((y - (slope_RMA_weighted * x + intercept_RMA_weighted))^2)))/(n - 2)) / (n * sum(weights * ((x - x_mean)^2))) )
    error_intercept_RMA_weighted <-  error_slope_RMA_weighted * sqrt(sum(weights * (x)^2) / (sum(weights))) 
    slope_MA_weighted <- (sum(weight_factors * v_w^2) - sum(weight_factors * u_w^2) + sqrt((sum(weight_factors * v_w^2) - sum(weight_factors * u_w^2))^2 + 4 * (sum(weight_factors * u_w * v_w)^2))) /
      (2 * sum(weight_factors * u_w * v_w))
    intercept_MA_weighted <- y_mean_weighted - slope_MA_weighted * x_mean_weighted
    error_slope_MA_weighted <- sqrt( ((n * sum(weights * ((y - (slope_MA_weighted * x + intercept_MA_weighted))^2)))/(n - 2)) / (n * sum(weights * ((x - x_mean)^2))) )
    error_intercept_MA_weighted <-  error_slope_MA_weighted * sqrt(sum(weights * (x)^2) / (sum(weights)))
    pearson_r = cor(x, y, method = "pearson")
    return(list(slopes = list(MA = slope_MA_normal, MA_weighted = slope_MA_weighted, RMA = slope_RMA_normal, RMA_weighted = slope_RMA_weighted), 
                errors_slopes = list(MA = error_slope_MA_normal, MA_weighted = error_slope_MA_weighted, RMA = error_slope_RMA_normal, RMA_weighted = error_slope_RMA_weighted), 
                intercepts = list(MA = intercept_MA_normal, MA_weighted = intercept_MA_weighted, RMA = intercept_RMA_normal, RMA_weighted = intercept_RMA_weighted),
                errors_intercepts = list(MA = error_intercept_MA_normal, MA_weighted = error_intercept_MA_weighted, RMA = error_intercept_RMA_normal, RMA_weighted = error_intercept_RMA_weighted),
                pearson_r = pearson_r))
  } else if (is.null(weights)) {
    n <- length(x)
    x_mean <- mean(x)
    y_mean <- mean(y)
    sum_x <- sum(x)
    sum_y <- sum(y)
    sum_x2 <- sum(x^2)
    sum_y2 <- sum(y^2)
    sum_u2 <- sum((x - x_mean)^2)
    sum_v2 <- sum((y - y_mean)^2)
    sum_uv <- sum((x - x_mean) * (y - y_mean))
    slope_RMA_normal <- sqrt((sum_y2 - sum_y^2/n) / (sum_x2 - sum_x^2/n))
    intercept_RMA_normal <- y_mean - slope_RMA_normal * x_mean
    error_slope_RMA_normal <- sqrt( ((sum((y - (slope_RMA_normal * x + intercept_RMA_normal))^2))/(n - 2)) / sum((x - x_mean)^2) )
    error_intercept_RMA_normal <- error_slope_RMA_normal * sqrt(sum_x2 / n)
    slope_MA_normal <- (sum_v2 - sum_u2 + sqrt((sum_v2 - sum_u2)^2 + 4 * ((sum_uv)^2))) / (2 * sum_uv)
    intercept_MA_normal <- y_mean - slope_MA_normal * x_mean
    error_slope_MA_normal <- sqrt( ((sum((y - (slope_MA_normal * x + intercept_MA_normal))^2))/(n - 2)) / sum((x - x_mean)^2) )
    error_intercept_MA_normal <- error_slope_MA_normal * sqrt(sum_x2 / n)
    pearson_r = cor(x, y, method = "pearson")
    return(list(slopes = list(MA = slope_MA_normal, RMA = slope_RMA_normal), 
                errors_slopes = list(MA = error_slope_MA_normal, RMA = error_slope_RMA_normal),
                intercepts = list(MA = intercept_MA_normal, RMA = intercept_RMA_normal),
                errors_intercepts = list(MA = error_intercept_MA_normal, RMA = error_intercept_RMA_normal),
                pearson_r = pearson_r))
  }
}

# Function to create one or multiple plots (one plot per facet) for the calculated regression
RMA_plot <- function(data, regression, x, y, facet_var = NULL) {
  p <- data %>%
    ggplot() +
    geom_point(aes(x = {{x}}, y = {{y}}), size = 0.8) +
    geom_abline(slope = regression$slopes$RMA_weighted, intercept = regression$intercepts$RMA_weighted, color = "darkblue") +
    geom_abline(slope = regression$slopes$RMA, intercept = regression$intercepts$RMA, color = "blue") +
    geom_abline(slope = regression$slopes$MA_weighted, intercept = regression$intercepts$MA_weighted, color = "darkred") +
    geom_abline(slope = regression$slopes$MA, intercept = regression$intercepts$MA, color = "red") +
    annotate("text", x = -13.3, y = 0,
             label = paste("PWRMA: y = ",
                           sprintf("%.2f", regression$slopes$RMA_weighted), " (", sprintf("%.2f", regression$errors_slopes$RMA_weighted), ") x + ",
                           sprintf("%.2f", regression$intercepts$RMA_weighted), " (", sprintf("%.2f", regression$errors_intercepts$RMA_weighted), ")",
                           sep = ""),
             size = 2, color = "darkblue") + 
    annotate("text", x = -13.3, y = -8,
             label = paste("RMA: y = ",
                           sprintf("%.2f", regression$slopes$RMA), " (", sprintf("%.2f", regression$errors_slopes$RMA), ") x + ",
                           sprintf("%.2f", regression$intercepts$RMA), " (", sprintf("%.2f", regression$errors_intercepts$RMA), ")",
                           sep = ""),
             size = 2, color = "blue") + 
    annotate("text", x = -13.3, y = -16,
             label = paste("PWMA: y = ",
                           sprintf("%.2f", regression$slopes$MA_weighted), " (", sprintf("%.2f", regression$errors_slopes$MA_weighted), ") x + ",
                           sprintf("%.2f", regression$intercepts$MA_weighted), " (", sprintf("%.2f", regression$errors_intercepts$MA_weighted), ")",
                           sep = ""),
             size = 2, color = "darkred") + 
    annotate("text", x = -13.3, y = -24,
             label = paste("MA: y = ",
                           sprintf("%.2f", regression$slopes$MA), " (", sprintf("%.2f", regression$errors_slopes$MA), ") x + ",
                           sprintf("%.2f", regression$intercepts$MA), " (", sprintf("%.2f", regression$errors_intercepts$MA), ")",
                           sep = ""),
             size = 2, color = "red") +
    annotate("text", x = -13.3, y = -32,
             label = paste("Pearson r:",
                           sprintf("%.4f", regression$pearson_r)),
             size = 2, color = "black") +
    labs(x = expression("\u03B4"^"18"*"O [\u2030]"),
         y = expression("\u03B4"^"2"*"H [\u2030]")) +
    scale_x_continuous(labels = ~sub("-", "\u2212", .x), limits = c(-20, 0)) +
    scale_y_continuous(labels = ~sub("-", "\u2212", .x), limits = c(-150, 0)) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 7, colour = "black"),
          axis.title.x = element_text(size = 7, colour = "black"),
          axis.text.y = element_text(size = 7, colour = "black"),
          axis.title.y = element_text(size = 7, colour = "black"),
          strip.text = element_text(size = 7, colour = "black"),
          legend.position = "none")
  return(p)
  p
}

RMA_plot_facet <- function(data, x, y, facet_var, regression) {
  RMA_plot(data, regression, {{x}}, {{y}}, {{facet_var}}) + facet_wrap(enquo(facet_var))
}

