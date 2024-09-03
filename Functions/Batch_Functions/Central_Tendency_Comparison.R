
# Central_Tendency_Comparison.R

# ONLY POSSIBLE IF TRUTH EXISTS

# truth_ppp: A spatstat ppp object, where each point is a true point corresponding to one pattern in a list of ppp objects
# centers_df: A dataframe of calculated centers from the list of ppp objects

compare_centers <- function(truth_ppp, centers_df) {
  
  # Distance calculations
  center_to_true_list <- list()
  loop_indices <- seq(from = 1, to = ncol(centers_df), by = 2)
  i <- 0
  for (column in loop_indices) {
    i <- i + 1
    center_to_true <- c()
    for (index in 1:nrow(centers_df)) {
      if (any(!is.na(truth_ppp[[index]]))) {
        x1y1 <- data.frame(X = centers_df[,column][index], Y = centers_df[,column + 1][index])
        x2y2 <- data.frame(X = truth_ppp[[index]]$x, Y = truth_ppp[[index]]$y)
        center_to_true <- c(center_to_true, sqrt(((x2y2$X - x1y1$X)^2 + (x2y2$Y - x1y1$Y)^2)))
      } else {
        center_to_true <- c(center_to_true, NA)
      }
    }
    center_to_true_list[[i]] <- center_to_true
  }
  distances_dataframe_withNA <- data.frame(matrix(NA, nrow = nrow(centers_df), ncol = (ncol(centers_df) / 2)))
  for (index in 1:length(center_to_true_list)) {
    distances_dataframe_withNA[,index] <- center_to_true_list[[index]]
  }
  store_colnames <- c()
  if (!is.na(which(colnames(centers_df) %in% "Mdx"))) {
    store_colnames <- c(store_colnames, "Distance_MD")
  }
  if (!is.na(which(colnames(centers_df) %in% "Mx"))) {
    store_colnames <- c(store_colnames, "Distance_MN")
  }
  if (!is.na(which(colnames(centers_df) %in% "Mwx"))) {
    store_colnames <- c(store_colnames, "Distance_WMN")
  }
  if (all(!is.na(which(grepl("C", colnames(centers_df)))))) {
    for (val in 1:length(which(grepl("C", colnames(centers_df)))[which(which(grepl("C", colnames(centers_df))) %% 2 == 1)])) {
      store_colnames <- c(store_colnames, paste0("Distance_WMNC_", val))
    }
  }
  colnames(distances_dataframe_withNA) <- store_colnames
  
  # Angle calculations
  alpha_list <- list()
  loop_indices <- seq(from = 1, to = ncol(centers_df), by = 2)
  i <- 0
  for (column in loop_indices) {
    i <- i + 1
    stored_alpha <- c()
    for (index in 1:nrow(centers_df)) {
      if (any(!is.na(truth_ppp[[index]]))) {
        y1 <- centers_df[,column + 1][index]
        yt <- truth_ppp[[index]]$y
        x1 <- centers_df[,column][index]
        xt <- truth_ppp[[index]]$x
        alpha <- atan2(yt - y1, xt - x1)
        if (alpha < 0) {
          alpha <- alpha + (2*pi)
        }
        stored_alpha <- c(stored_alpha, rad2deg(alpha))
      } else {
        stored_alpha <- c(stored_alpha, NA)
      }
    }
    alpha_list[[i]] <- stored_alpha
  }
  angles_dataframe_withNA <- data.frame(matrix(NA, nrow = nrow(centers_df), ncol = (ncol(centers_df) / 2)))
  for (index in 1:length(alpha_list)) {
    angles_dataframe_withNA[,index] <- alpha_list[[index]]
  }
  new_strings <- unlist(strsplit(colnames(distances_dataframe_withNA), split='Distance_', fixed=TRUE))
  new_strings <- new_strings[nzchar(new_strings)]
  colnames(angles_dataframe_withNA) <- paste0("Angle_", new_strings)
  
  # Dataframe compilation
  distances_dataframe <- na.omit(distances_dataframe_withNA)
  angles_dataframe <- na.omit(angles_dataframe_withNA)
  
  # Comparison of measures
  all_comparisons <- t(combn(1:ncol(distances_dataframe), 2))
  
  # y-axis is further from true P% of the time than x-axis
  comp_matrix <- matrix(NA, nrow = length(distances_dataframe), ncol = length(distances_dataframe))
  for (rows in 1:nrow(all_comparisons)) {
    comp_matrix[all_comparisons[rows,][1], all_comparisons[rows,][2]] <- sum(distances_dataframe[,all_comparisons[rows,][1]] > distances_dataframe[,all_comparisons[rows,][2]]) / nrow(distances_dataframe)
    comp_matrix[all_comparisons[rows,][2], all_comparisons[rows,][1]] <- sum(distances_dataframe[,all_comparisons[rows,][2]] > distances_dataframe[,all_comparisons[rows,][1]]) / nrow(distances_dataframe)
  }
  comp_matrix[is.na(comp_matrix)] <- 0
  colnames(comp_matrix) <- new_strings
  rownames(comp_matrix) <- new_strings
  
  # When y-axis is further from true than x-axis, it is, on average, P% further
  comp_matrix_2 <- matrix(NA, nrow = length(distances_dataframe), ncol = length(distances_dataframe))
  for (rows in 1:nrow(all_comparisons)) {
    comp_matrix_2[all_comparisons[rows,][1], all_comparisons[rows,][2]] <- (sum(distances_dataframe[,all_comparisons[rows,][1]][distances_dataframe[,all_comparisons[rows,][1]] > distances_dataframe[,all_comparisons[rows,][2]]]) / sum(distances_dataframe[,all_comparisons[rows,][2]][distances_dataframe[,all_comparisons[rows,][1]] > distances_dataframe[,all_comparisons[rows,][2]]])) - 1
    comp_matrix_2[all_comparisons[rows,][2], all_comparisons[rows,][1]] <- (sum(distances_dataframe[,all_comparisons[rows,][2]][distances_dataframe[,all_comparisons[rows,][2]] > distances_dataframe[,all_comparisons[rows,][1]]]) / sum(distances_dataframe[,all_comparisons[rows,][1]][distances_dataframe[,all_comparisons[rows,][2]] > distances_dataframe[,all_comparisons[rows,][1]]])) - 1
  }
  comp_matrix_2[is.na(comp_matrix_2)] <- 0
  colnames(comp_matrix_2) <- new_strings
  rownames(comp_matrix_2) <- new_strings
  
  comparison_list <- list(comp_matrix, comp_matrix_2, distances_dataframe_withNA, angles_dataframe_withNA)
  names(comparison_list) <- c("distance_pct", "magnitude_pct", "distances", "angles")
  
  return(comparison_list)
  
} 
