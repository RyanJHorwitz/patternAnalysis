
# Central_Tendency_Metrics.R

# -- REQUIRED
# measures: Vector of selected measures of central tendency. Options: "all" or c("mean", "weighted_mean", "constrained_weighted_mean")
# distance_ppp: list of spatstat ppp objects with equidistant projection
# storage_path: Where to save outputs of function.  If null, new folder created in cwd
# overwrite: Boolean, whether to overwrite saved outputs of function from prior run

# -- When measures are "all" or "weighted_mean_center"
# area_ppp: list of spatstat ppp objects with equal-area projection
# neighborhood_size: Number of nearest neighbors used in weighting calculation. Weights are calculated for each point in an input pattern, such that weight equals the neighborhood size (N) divided by the area of the minimum bounding rectangle around point (P) and its N nearest neighbors.

# -- When measures are "all" or "constrained_weighted_mean"
# area_ppp: list of spatstat ppp objects with equal-area projection
# area_crs: An equal area projection
# neighborhood_size: Number of nearest neighbors used in weighting calculation. Weights are calculated for each point, such that weight equals the neighborhood size (N) divided by the area of the minimum bounding rectangle around point (P) and its N nearest neighbors.
# constraint_paths: Filepath (or vector of filepaths) to known spatial range of species of interest

calculate_centers <- function(measures, area_crs, area_ppp, distance_ppp, neighborhood_size = 50, constraint_paths, storage_path = NULL, overwrite = FALSE) {
  
  if (is.null(storage_path)) {
    storage_path <- "./Centers_Function_Checkpoints"
    dir.create(storage_path)
  }
  
  if (measures == "all") {
    measures <- c("median", "mean", "weighted_mean", "constrained_weighted_mean")
  }
  
  if (all(grepl("constrained_weighted_mean", measures) & !grepl("weighted_mean", measures))) {
    measures <- c(measures, "weighted_mean")
  }
  
  # MEDIAN CENTER
  if ("median" %in% measures) {
    
    # Minimize the Euclidean Distance to all points for each elephant point pattern 
    
    all_ppp <- distance_ppp
    
    median_x <- c()
    median_y <- c()
    for (ppp_index in 1:length(all_ppp)) {
      x <- coords(all_ppp[[ppp_index]])[,1]
      y <- coords(all_ppp[[ppp_index]])[,2]
      crit <- 0.0001
      x0 <- mean(x)
      y0 <- mean(y)
      dx <- Inf
      dy <- Inf
      iteration <- 0
      while (abs(dx) > crit | abs(dy) > crit) {
        xd <- x - x0
        yd <- y - y0
        d <- sqrt(xd*xd + yd*yd)
        w <- as.numeric(1/d)
        w <- w / sum(w)
        x1 <- w * x
        x1 <- sum(x1)
        y1 <- w * y
        y1 <- sum(y1)
        dx <- x1 - x0
        dy <- y1 - y0
        iteration <- iteration + 1
        x0 <- x1
        y0 <- y1
      }
      output <- c(x1, y1)
      median_x <- c(median_x, output[1])
      median_y <- c(median_y, output[2])
    }
  }
  
  # MEAN CENTER
  if ("mean" %in% measures) {
    
    # Take the average of all x/y coordinates for each elephant point pattern
    
    all_ppp <- distance_ppp
    
    mean_x <- c()
    mean_y <- c()
    for (ppp_index in 1:length(all_ppp)) {
      mean_x <- c(mean_x, mean(coords(all_ppp[[ppp_index]])[,1])) # Mean x-coordinate
      mean_y <- c(mean_y, mean(coords(all_ppp[[ppp_index]])[,2])) # Mean y-coordinate
    }
  }
  
  # WEIGHTED MEAN CENTER
  if ("weighted_mean" %in% measures) {
    
    ### For each elephant sample, compute the matrix of distances between all points associated with the sample
    ### For each row in the matrix (where each row represents a single point, "p"), identify the [neighborhood size] points closest to "p"
    ### Extract the identified closest points, including the row sample and draw the minimum bounding geometry needed to enclose the points
    ### Divide [neighborhood size + 1] points by the area of the minimum bounding geometry
    ### Store in list of dataframes where each dataframe refers to all point weights of a single elephant sample
    
    all_ppp <- area_ppp # Using equal area projection for density calculation
    all_ppp_equidistant <- distance_ppp
    smooth <- 10e-15 # PREVENTS DIVIDE BY ZERO ERROR IN CASE ALL POINTS IN NEIGHBORHOOD ARE PREDICTED IN SAME LOCATION
    
    if (!file.exists(paste0(storage_path, "stored_top_all_ppp.RDS")) | !file.exists(paste0(storage_path, "stored_pattern_all_ppp.RDS")) | overwrite == TRUE) {
      stored_top_all_ppp <- list()
      stored_pattern_all_ppp <- list()
      for (ppp_index in 1:length(all_ppp)) {
        dataset <- pairdist(all_ppp_equidistant[[ppp_index]])
        stored_top <- list()
        pattern <- data.frame()
        sample_size <- npoints(all_ppp[[ppp_index]])
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
          extraction <- coords(all_ppp[[ppp_index]])[top,]
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
          print(paste0(row, " / ", nrow(dataset), " || ", ppp_index, " / ", length(all_ppp)))
        }
        stored_top_all_ppp[[ppp_index]] <- stored_top
        stored_pattern_all_ppp[[ppp_index]] <- pattern
      }
      saveRDS(stored_pattern_all_ppp, paste0(storage_path, "stored_pattern_all_ppp.RDS"))
      saveRDS(stored_top_all_ppp, paste0(storage_path, "stored_top_all_ppp.RDS"))
    } else {
      stored_pattern_all_ppp <- readRDS(paste0(storage_path, "stored_pattern_all_ppp.RDS"))
      stored_top_all_ppp <- readRDS(paste0(storage_path, "stored_top_all_ppp.RDS"))
    }
    
    # Switch to equidistant projection as central tendency measures are distance-based
    all_ppp <- distance_ppp # Note change in coordinate system to equidistant
    
    # Store each elephant sample in a list of dataframes, where each dataframe includes the coordinates of all 
    # points in single elephant sample's point pattern and their associated weights
    # Calculate the median, mean, and weighted mean for each sample
    # Store all calculated centers in a single, combined dataframe
    
    if (!file.exists(paste0(storage_path, "stored_full_dataset_all_ppp.RDS")) | overwrite == TRUE) {
      stored_full_dataset_all_ppp <- list()
      w_mean_x <- c()
      w_mean_y <- c()
      for (ppp_index in 1:length(all_ppp)) {
        x <- all_ppp[[ppp_index]]$x
        y <- all_ppp[[ppp_index]]$y
        rect <- stored_pattern_all_ppp[[ppp_index]]$rect_dens
        circ <- stored_pattern_all_ppp[[ppp_index]]$circ_dens
        full_dataset <- data.frame(x, y, rect, circ)
        stored_full_dataset_all_ppp[[ppp_index]] <- full_dataset
        w_mean_x <- c(w_mean_x, sum(full_dataset$x * full_dataset$rect) / sum(full_dataset$rect)) # Weighted mean x-coordinate
        w_mean_y <- c(w_mean_y, sum(full_dataset$y * full_dataset$rect) / sum(full_dataset$rect)) # Weighted mean y-coordinate
      }
      stored_w_mean_coordinates <- data.frame(w_mean_x = w_mean_x, w_mean_y = w_mean_y)
      saveRDS(stored_w_mean_coordinates, paste0(storage_path, "stored_w_mean_coordinates.RDS"))
      saveRDS(stored_full_dataset_all_ppp, paste0(storage_path, "stored_full_dataset_all_ppp.RDS"))
    } else {
      stored_w_mean_coordinates <- readRDS(paste0(storage_path, "stored_w_mean_coordinates.RDS"))
      stored_full_dataset_all_ppp <- readRDS(paste0(storage_path, "stored_full_dataset_all_ppp.RDS"))
    }
  }
  
  # CONSTRAINED WEIGHTED MEAN CENTER
  if ("constrained_weighted_mean" %in% measures) {
    
    constraints_to_test <- list()
    for (path_index in 1:length(constraint_paths)) {
      constraint <- st_geometry(read_sf(dsn = constraint_paths[path_index])) # Read in data
      constraint_projected <- st_transform(constraint, crs = area_crs) # Project coordinates
      constraints_to_test[[path_index]] <- constraint_projected
    }
    
    # For all constraints, for each elephant point pattern in the dataset...
    # -- Find the points that fall outside of the constraint
    # -- Down weight the points to 0
    # -- Recalculate the weighted mean (now called the constrained weighted mean)
    # -- Store the new constrained weighted means in a list of dataframes where... 
    # ----- Each list index corresponds to a constraint
    # ----- Each dataframe contains the constrained weighted mean coordinates for all elephant samples in the dataset
    # -- Store the point reweights in a list of lists of dataframes where...
    # ----- Each outer list index corresponds to a constraint
    # ----- Each inner list index corresponds to an elephant sample point pattern
    # ----- Each dataframe corresponds to the reweights and coordinates of all points in a single elephant sample
    
    if (!file.exists(paste0(storage_path, "constrained_weighted_list.RDS")) | !file.exists(paste0(storage_path, "constrained_reweights_list.RDS"))) {
      constrained_weighted_list <- list()
      constrained_reweights_list <- list()
      for (sf_index in 1:length(constraints_to_test)) {
        constrained_weighted <- data.frame()
        reweights_list <- list()
        for (elephants in 1:length(stored_full_dataset_all_ppp)) {
          points_sf <- st_as_sf(stored_full_dataset_all_ppp[[elephants]], coords = c("x", "y"), crs = area_crs)
          st_agr(points_sf) = "constant" # Attribute variables are assumed to be spatially constant throughout all geometries
          point_intersection <- st_intersection(points_sf, constraints_to_test[[sf_index]])
          intersection_df <- data.frame(cbind(dplyr::select(as.data.frame(point_intersection), -geometry), st_coordinates(point_intersection)))
          cw_mean_x <- sum(intersection_df$X * intersection_df$rect) / sum(intersection_df$rect)
          cw_mean_y <- sum(intersection_df$Y * intersection_df$rect) / sum(intersection_df$rect)
          constrained_weighted <- rbind(constrained_weighted, data.frame(x = cw_mean_x, y = cw_mean_y))
          
          point_difference <- st_difference(points_sf, constraints_to_test[[sf_index]])
          difference_df <- data.frame(cbind(dplyr::select(as.data.frame(point_difference), -geometry), st_coordinates(point_difference)))
          difference_df$rect <- rep(0, nrow(difference_df))
          difference_df$circ <- rep(0, nrow(difference_df))
          
          reweights <- rbind(intersection_df, difference_df)
          reweights$index <- as.numeric(row.names(reweights))
          reweights <- reweights[order(reweights$index), ]
          reweights_sorted <- data.frame(rect = reweights$rect, circ = reweights$circ, x = reweights$X, y = reweights$Y)
          reweights_list[[elephants]] <- reweights_sorted
          
          print(paste0(elephants, " / ", length(stored_full_dataset_all_ppp), " || ", sf_index, " / ", length(constraints_to_test)))
        }
        constrained_weighted_list[[sf_index]] <- constrained_weighted
        constrained_reweights_list[[sf_index]] <- reweights_list
      }
      saveRDS(constrained_weighted_list, paste0(storage_path, "constrained_weighted_list.RDS"))
      saveRDS(constrained_reweights_list, paste0(storage_path, "constrained_reweights_list.RDS"))
    } else {
      constrained_weighted_list <- readRDS(paste0(storage_path, "constrained_weighted_list.RDS"))
      constrained_reweights_list <- readRDS(paste0(storage_path, "constrained_reweights_list.RDS"))
    }
  }
  
  if (!file.exists(paste0(storage_path, "stored_centers_all_ppp.RDS")) | overwrite == TRUE) {
    stored_centers_all_ppp <- data.frame(rep(NA, length(all_ppp)))
    for (metrics in measures) {
      if ("median" %in% metrics) {
        median_df <- data.frame(Mdx = median_x, Mdy = median_y)
        stored_centers_all_ppp <- data.frame(stored_centers_all_ppp, median_df)
      }
      if ("mean" %in% metrics) {
        mean_df <- data.frame(Mx = mean_x, My = mean_y)
        stored_centers_all_ppp <- data.frame(stored_centers_all_ppp, mean_df)
      }
      if ("weighted_mean" %in% metrics) {
        weighted_mean_df <- data.frame(Mwx = stored_w_mean_coordinates$w_mean_x, Mwy = stored_w_mean_coordinates$w_mean_y)
        stored_centers_all_ppp <- data.frame(stored_centers_all_ppp, weighted_mean_df)
      }
      if ("constrained_weighted_mean" %in% metrics) {
        constrained_weighted_mean_df <- data.frame(rep(NA, length(all_ppp)))
        for (index in 1:length(constrained_weighted_list)) {
          constrained_weighted_mean_df <- data.frame(constrained_weighted_mean_df, constrained_weighted_list[[index]])
        }
        constrained_weighted_mean_df <- constrained_weighted_mean_df[-1]
        prep_cols <- rep(c(paste0("Mwx_C"), paste0("Mwy_C")), 2)
        prep_inds <- sort(rep(1:(ncol(constrained_weighted_mean_df) / 2), 2))
        colnames(constrained_weighted_mean_df) <- paste0(prep_cols, "_", prep_inds)
        stored_centers_all_ppp <- data.frame(stored_centers_all_ppp, constrained_weighted_mean_df)
      }
    }
    stored_centers_all_ppp <- stored_centers_all_ppp[,-1]
    saveRDS(stored_centers_all_ppp, paste0(storage_path, "stored_centers_all_ppp.RDS"))
  } else {
    stored_centers_all_ppp <- readRDS(paste0(storage_path, "stored_centers_all_ppp.RDS"))
  }
  
  if (all(!grepl("weighted_mean", measures) & !grepl("constrained_weighted_mean", measures))) {
    centers_outputs <- list(measures, stored_centers_all_ppp)
    names(centers_outputs) <- c("measures", "centers")
    return(centers_outputs)
  } else if (all(grepl("weighted_mean", measures) & !grepl("constrained_weighted_mean", measures))) {
    centers_outputs <- list(measures, stored_centers_all_ppp, stored_full_dataset_all_ppp)
    names(centers_outputs) <- c("measures", "centers", "weights")
    return(centers_outputs)
  } else {
    centers_outputs <- list(measures, stored_centers_all_ppp, stored_full_dataset_all_ppp, constrained_reweights_list)
    names(centers_outputs) <- c("measures", "centers", "weights", "constrained_weights")
    return(centers_outputs)
  }


}