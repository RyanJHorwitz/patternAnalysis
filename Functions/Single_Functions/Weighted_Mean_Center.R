
# Weighted_Mean_Center.R

# Source: https://pro.arcgis.com/en/pro-app/latest/tool-reference/spatial-statistics/h-how-mean-center-spatial-statistics-works.htm

# area_ppp: A spatstat ppp object with an equalarea projection
# distance_ppp: A spatstat ppp object with an equidistant projection
# neighborhood_size: Number of nearest neighbors used in weighting calculation. Weights are calculated for each point in an input pattern, such that weight equals the neighborhood size (N) divided by the area of the minimum bounding rectangle around point (P) and its N nearest neighbors.

weighted_mean_center <- function (area_ppp, distance_ppp, neighborhood_size = 50) {
  
  ### Compute the matrix of distances between all points associated with the sample
  ### For each row in the matrix (where each row represents a single point, "p"), identify the [neighborhood size] points closest to "p"
  ### Extract the identified closest points, including the row sample and draw the minimum bounding geometry needed to enclose the points
  ### Divide [neighborhood size + 1] points by the area of the minimum bounding geometry
  
  smooth <- 10e-15 # PREVENTS DIVIDE BY ZERO ERROR IN CASE ALL POINTS IN NEIGHBORHOOD ARE PREDICTED IN SAME LOCATION

  stored_top_all_ppp <- list()
  stored_pattern_all_ppp <- list()

  dataset <- pairdist(distance_ppp)
  stored_top <- list()
  pattern <- data.frame()
  sample_size <- npoints(area_ppp)
  for (row in 1:nrow(dataset)) {
    point <- dataset[row,]
    names(point) <- 1:sample_size
    sorted_points <- sort(point)
    if (neighborhood_size >= 1) {
      top <- as.numeric(names(sorted_points)[1:(neighborhood_size + 1)]) # Includes itself, hence add 1
    } else {
      top <- as.numeric(names(sorted_points)[1:(as.integer(length(sorted_points) * neighborhood_size) + 1)]) # Includes itself, hence add 1
    }
    stored_top[[row]] <- top
    extraction <- coords(area_ppp)[top,]
    min_bounding_rect <- boundingbox(extraction)
    if (area.owin(as.owin(min_bounding_rect)) == 0) { # IN CASE ALL NEIGHBORHOOD POINTS ARE SAME LOCATION
      min_bounding_circ <- smooth
      if (neighborhood_size >= 1) {
        rect_dens <- (neighborhood_size + 1) / smooth
        circ_dens <- (neighborhood_size + 1) / smooth
      } else {
        rect_dens <- (as.integer(length(sorted_points) * neighborhood_size) + 1) / smooth
        circ_dens <- (as.integer(length(sorted_points) * neighborhood_size) + 1) / smooth
      }
    } else {
      min_bounding_circ <- boundingcircle(as.ppp(extraction, W = min_bounding_rect))
      if (neighborhood_size >= 1) {
        rect_dens <- ((neighborhood_size + 1) / area.owin(as.owin(min_bounding_rect))) + smooth
        circ_dens <- ((neighborhood_size + 1) / area.owin(as.owin(min_bounding_circ))) + smooth
      } else {
        rect_dens <- ((as.integer(length(sorted_points) * neighborhood_size) + 1) / area.owin(as.owin(min_bounding_rect))) + smooth
        circ_dens <- ((as.integer(length(sorted_points) * neighborhood_size) + 1) / area.owin(as.owin(min_bounding_circ))) + smooth
      }
    }
    pattern <- data.frame(rbind(pattern, (data.frame(rect_dens, circ_dens))))
  }
  
  x <- distance_ppp$x
  y <- distance_ppp$y
  rect <- pattern$rect_dens
  circ <- pattern$circ_dens
  full_dataset <- data.frame(x, y, rect, circ)
  stored_full_dataset_area_ppp <- full_dataset
  w_mean_x <- sum(full_dataset$x * full_dataset$rect) / sum(full_dataset$rect) # Weighted mean x-coordinate
  w_mean_y <- sum(full_dataset$y * full_dataset$rect) / sum(full_dataset$rect) # Weighted mean y-coordinate
      
  output <- list(coords = data.frame(w_mean_x = w_mean_x, w_mean_y = w_mean_y), area_values = data.frame(rect = full_dataset$rect, circ = full_dataset$circ))
  
  return (output)
}

