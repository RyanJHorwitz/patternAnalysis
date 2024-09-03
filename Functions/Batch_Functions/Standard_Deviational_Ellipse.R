
# General_Explanation: https://pro.arcgis.com/en/pro-app/latest/tool-reference/spatial-statistics/h-how-directional-distribution-standard-deviationa.htm
# Weighted_Variance: chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/http://seismo.berkeley.edu/~kirchner/Toolkits/Toolkit_12.pdf
# Weighted_Covariance: https://www.itl.nist.gov/div898/software/dataplot/refman2/auxillar/weigcorr.htm

# pattern_dataset: A list of spatstat ppp objects
# coordinates_colnames: Vector of column names (as string) of the X and Y coordinates of the ppp objects, respectively
# chosen_center: A spatstat ppp object containing the XY coordinates representing the center of the point pattern. For instance, mean center, median center, etc.
# coordinates_colnames: Vector of column names (as string) of the X and Y coordinates of the chosen_center ppp object, respectively
# weights_colname: The column name of the weights column (if applicable)
# adjustments: number of (weighted) standard distances/deviations from the (weighted) mean

standard_dev_ellipse <- function(pattern_dataset, coordinates_colnames, chosen_center, center_colnames, weights_colname = NA, adjustments = 1) {
  
  if (!is.na(weights_colname)) {
    # Calculation of the weighted variance and covariance
    var_x_ppp <- c()
    var_y_ppp <- c()
    covariance <- c()
    for (index in 1:length(pattern_dataset)) {
      n <- nrow(pattern_dataset[[index]]) # Number of points in pattern
      xm <- chosen_center[index,][,colnames(chosen_center[index,]) == center_colnames[1]] # x-coordinate of center
      ym <- chosen_center[index,][,colnames(chosen_center[index,]) == center_colnames[2]] # y-coordinate of center
      summed_top_x <- 0
      summed_top_y <- 0
      cov_summed_top <- 0
      summed_bottom <- 0
      for (point in 1:nrow(pattern_dataset[[index]])) {
        weight <- pattern_dataset[[index]][,colnames(pattern_dataset[[index]]) == weights_colname][point] # weight of point
        xi <- pattern_dataset[[index]][,colnames(pattern_dataset[[index]]) == coordinates_colnames[1]][point] # x-coordinate of point
        summed_top_x <- summed_top_x + weight*(xi - xm)^2
        yi <- pattern_dataset[[index]][,colnames(pattern_dataset[[index]]) == coordinates_colnames[2]][point] # y-coordinate of point
        summed_top_y <- summed_top_y + weight*(yi - ym)^2
        cov_summed_top <- cov_summed_top + (weight*(xi - xm) * (yi-ym))
        summed_bottom <- summed_bottom + weight
      }
      var_x_ppp <- c(var_x_ppp, (summed_top_x / ((summed_bottom * (n - 1)) / n)))
      var_y_ppp <- c(var_y_ppp, (summed_top_y / ((summed_bottom * (n - 1)) / n)))
      covariance <- c(covariance, (cov_summed_top / summed_bottom))
    }
    
  } else {
    
    # Calculation of the weighted variance and covariance
    var_x_ppp <- c()
    var_y_ppp <- c()
    covariance <- c()
    for (index in 1:length(pattern_dataset)) {
      n <- nrow(pattern_dataset[[index]]) # Number of points in pattern
      xm <- chosen_center[index,][,colnames(chosen_center[index,]) == center_colnames[1]] # x-coordinate of center
      ym <- chosen_center[index,][,colnames(chosen_center[index,]) == center_colnames[2]] # y-coordinate of center
      summed_top_x <- 0
      summed_top_y <- 0
      cov_summed_top <- 0
      summed_bottom <- 0
      for (point in 1:nrow(pattern_dataset[[index]])) {
        xi <- pattern_dataset[[index]][,colnames(pattern_dataset[[index]]) == coordinates_colnames[1]][point] # x-coordinate of point
        summed_top_x <- summed_top_x + (xi - xm)^2
        yi <- pattern_dataset[[index]][,colnames(pattern_dataset[[index]]) == coordinates_colnames[2]][point] # y-coordinate of point
        summed_top_y <- summed_top_y + (yi - ym)^2
        cov_summed_top <- cov_summed_top + ((xi - xm) * (yi - ym))
      }
      var_x_ppp <- c(var_x_ppp, ((summed_top_x / (n - 1))))
      var_y_ppp <- c(var_y_ppp, ((summed_top_y / (n - 1))))
      covariance <- c(covariance, ((cov_summed_top / n)))
    }
  }
  
  stored_major <- c()
  stored_minor <- c()
  stored_theta <- c()
  stored_ellipses <- list()
  for (index in 1:length(pattern_dataset)) {
    covariance_matrix <- as.matrix(rbind(c(var_x_ppp[index], covariance[index]), c(covariance[index], var_y_ppp[index])))
    a <- covariance_matrix[1, 1]
    b <- covariance_matrix[1, 2]
    c <- covariance_matrix[2, 2]
    eigens <- eigen(covariance_matrix) # Eigenvalues and Eigenvectors for covariance matrix
    angles <- seq(0, 2*pi, length.out=200) # Output 360 points from 0 - 2pi angles
    # Calculation of angle of rotation of ellipse (theta)
    if (b == 0 & a >= c) {
      theta <- 0
    } else if (b == 0 & a < c) {
      theta <- pi/2
    } else {
      theta <- atan2(eigens$values[1] - a, b)
    }
    stored_theta <- c(stored_theta, theta)
    # Calculation of ellipse coordinates as points
    ell <- eigens$vectors %*% t(cbind(adjustments*sqrt(eigens$values[1])*cos(angles), adjustments*sqrt(eigens$values[2])*sin(angles)))
    stored_ellipses[[index]] <- data.frame(x = chosen_center[index,][,colnames(chosen_center[index,]) == center_colnames[1]] + ell[1,], y = chosen_center[index,][,colnames(chosen_center[index,]) == center_colnames[2]] + ell[2,])
    # Calculation of major axis
    stored_major <- c(stored_major, adjustments*sqrt(eigens$values[1]))
    # Calculation of minor axis
    stored_minor <- c(stored_minor, adjustments*sqrt(eigens$values[2]))
    print(paste0(index, " / ", length(pattern_dataset)))
  }
  
  ellipse <- list(stored_ellipses, stored_major, stored_minor, stored_theta)
  names(ellipse) <- c("ellipses", "major", "minor", "theta")
  return(ellipse)
  
  return(covariance)
}