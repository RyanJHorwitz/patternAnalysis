
# Point_Densities.R

# input_sf: Input sf object with equal-area projection
# bounding_geometry: Set to "circle" or "rectangle"

global_density <- function(input_sf, bounding_geometry = "rectangle") {
  ppp_ea <- as.ppp(input_sf, W = as.owin(st_bbox(input_sf)))
  if (bounding_geometry == "rectangle") {
    min_bounding <- boundingbox(ppp_ea)
  }
  if (bounding_geometry == "circle") {
    min_bounding <- boundingcircle(ppp_ea)
  }
  if (area.owin(as.owin(min_bounding)) == 0) { # IN CASE ALL NEIGHBORHOOD POINTS ARE SAME LOCATION
    stop("All points overlap or bounding geometry failed to form around points.")
  } else {
    gde <- npoints(ppp_ea) / area.owin(as.owin(min_bounding))
  }
  return (gde)
}

# input_equidistant: Input sf object with equidistant projection
# input_equalarea: Input sf object with equal-area projection
# neighborhood_size: If less than 1, set to percentage of total # of points. Else, set to absolute number of points
# bounding_geometry: Set to "circle" or "rectangle"
# ordered_samples: Vector of ordered sample names to append to output densities

local_density <- function(input_equidistant, input_equalarea, neighborhood_size = 0.5, bounding_geometry = "rectangle", ordered_samples = NULL) {
  ppp_ed <- as.ppp(input_equidistant, W = as.owin(st_bbox(input_equidistant)))
  ppp_ea <- as.ppp(input_equalarea, W = as.owin(st_bbox(input_equalarea)))
  dataset <- pairdist(ppp_ed)
  stored_top <- list()
  pattern <- c()
  sample_size <- npoints(ppp_ed)
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
    extraction <- coords(ppp_ea)[top,]
  
    if (bounding_geometry == "rectangle") {
      min_bounding <- boundingbox(extraction)
    }
    if (bounding_geometry == "circle") {
      min_bounding <- boundingcircle(extraction)
    }
    if (area.owin(as.owin(min_bounding)) == 0) { # IN CASE ALL NEIGHBORHOOD POINTS ARE SAME LOCATION
      stop("All points overlap or bounding geometry failed to form around points.")
    } else {
      if (neighborhood_size >= 1) {
        dens <- ((neighborhood_size + 1) / area.owin(as.owin(min_bounding)))
      } else {
        dens <- ((as.integer(length(sorted_points) * neighborhood_size) + 1) / area.owin(as.owin(min_bounding)))
      }
    }
      pattern <- c(pattern, dens)
  }
  if (!is.null(ordered_samples)) {
    output <- data.frame(samples = ordered_samples, densities = pattern)
  } else {
    output <- data.frame(densities = pattern)
  }
  return(output)
}