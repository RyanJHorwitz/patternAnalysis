
# Constrained_Weighted_Mean_Center.R

# Source: https://pro.arcgis.com/en/pro-app/latest/tool-reference/spatial-statistics/h-how-mean-center-spatial-statistics-works.htm
# Constraint adaptation is custom

# area_ppp: A spatstat ppp object with an equalarea projection
# distance_ppp: A spatstat ppp object with an equidistant projection
# neighborhood_size: Number of nearest neighbors used in weighting calculation. Weights are calculated for each point in an input pattern, such that weight equals the neighborhood size (N) divided by the area of the minimum bounding rectangle around point (P) and its N nearest neighbors.
# constraint_paths: Filepath (or vector of filepaths) to known spatial range of species of interest

constrained_weighted_mean <- function(area_ppp, distance_ppp, area_crs, neighborhood_size = 50, constraint_paths) {
  
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

  rect <- pattern$rect_dens
  circ <- pattern$circ_dens
  wm <- data.frame(rect, circ)
  
  constraints_to_test <- list()
  for (path_index in 1:length(constraint_paths)) {
    constraint <- st_geometry(read_sf(dsn = constraint_paths[path_index])) # Read in data
    constraint_projected <- st_transform(constraint, crs = area_crs) # Project coordinates
    constraints_to_test[[path_index]] <- constraint_projected
  }
  
  # For all constraints...
  # -- Find the points that fall outside of the constraint
  # -- Down weight the points to 0
  # -- Recalculate the weighted mean (now called the constrained weighted mean)
  # -- Store the new constrained weighted means in a list of dataframes where... 
  # ----- Each list index corresponds to a constraint
  # ----- Each dataframe contains the constrained weighted mean coordinates for the elephant sample in the dataset
  # -- Store the point reweights in a list of dataframes where...
  # ----- Each list index corresponds to a constraint
  # ----- Each dataframe corresponds to the reweights and coordinates of all points in the elephant sample

  constrained_weighted_list <- list()
  constrained_reweights_list <- list()
  for (sf_index in 1:length(constraints_to_test)) {
    constrained_weighted <- data.frame()
    reweights_list <- list()
    
    points_sf <- st_as_sf(area_ppp, coords = c("x", "y"))[-1,][[2]]
    st_crs(points_sf) <- area_crs
    points_sf <- st_sf(points_sf)
    
    st_agr(points_sf) = "constant" # Attribute variables are assumed to be spatially constant throughout all geometries
    st_crs(points_sf) <- area_crs
    point_intersection <- st_intersection(points_sf, constraints_to_test[[sf_index]])
    intersection_df <- data.frame(cbind(dplyr::select(as.data.frame(point_intersection), -"points_sf"), st_coordinates(point_intersection)))
    cw_mean_x <- sum(intersection_df$X * wm$rect) / sum(wm$rect)
    cw_mean_y <- sum(intersection_df$Y * wm$rect) / sum(wm$rect)
    constrained_weighted <- data.frame(cw_mean_x = cw_mean_x, cw_mean_y = cw_mean_y)

    constrained_weighted_list[[sf_index]] <- constrained_weighted
    
  }
  
  output <- constrained_weighted_list
  
  return (output)
  
}
