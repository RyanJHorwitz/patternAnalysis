
# FILE: Pattern_Analysis.R

# ---------------- SET WORKING DIRECTORY ---------------- #

library(rstudioapi)
setwd(sub('/[^/]*$', '', getSourceEditorContext()$path))


# ------------------- IMPORT PACKAGES ------------------- #

library(sf)
library(PROJ)
library(stringr)
library(dplyr)
library(spatstat)
library(pracma) # rad2deg() function
library(paletteer) # map gradient
library(DescTools) # DrawEllipse() function


# --------------- IMPORT CUSTOM FUNCTIONS --------------- #

source("../Functions/Batch_Functions/Central_Tendency_Metrics.R")
source("../Functions/Batch_Functions/Central_Tendency_Comparison.R")
source("../Functions/Batch_Functions/Standard_Deviational_Ellipse.R")
source("../Functions/Batch_Functions/Standard_Distance.R")


# --- READ-IN SETTINGS AND DATASET FROM CONFIG.R FILE --- #

source("./Settings/config.R")


# --------- CALCULATE CENTRAL TENDENCY MEASURES -------- #

# Call "calculate_centers" function
all_centers <- calculate_centers(measures = "all",
                                 area_crs = equalarea,
                                 area_ppp = ppp_list_equalarea, 
                                 distance_ppp = ppp_list_equidistant, 
                                 neighborhood_size = neighborhood_size, 
                                 constraint_paths = constraint_paths,
                                 storage_path = variables_path,
                                 overwrite = FALSE
)

# Store output of function as individual variables
measures <- all_centers$measures
stored_full_dataset_all_ppp <- all_centers$weights
stored_centers_all_ppp <- all_centers$centers
constrained_reweights_list <- all_centers$constrained_weights

# If data has sample groupings (e.g. seizures), add in group information to centers list
if (run_type == "real" | run_type == "simulated") {
  stored_centers_all_ppp$groupID <- rep(basename(input_list)[which(stored_group_lengths > 0)], stored_group_lengths[stored_group_lengths > 0])
  split_centers <- split(stored_centers_all_ppp, f = factor(stored_centers_all_ppp$groupID, levels = unique(stored_centers_all_ppp$groupID)))
} else if (run_type == "drop") {
  stored_centers_all_ppp$groupID <- drop_zones
  split_centers <- split(stored_centers_all_ppp, f = factor(stored_centers_all_ppp$groupID, levels = unique(stored_centers_all_ppp$groupID)))
}


# ---------- COMPARE CENTRAL TENDENCY MEASURES --------- #

# Quantitative comparison (Only possible if truth.exists)
if (truth.exists) {
  if (run_type == "simulated" | run_type == "drop") {
    center_comparison <- compare_centers(truth_ppp = equidistant_truth_list_ppp, centers_df = stored_centers_all_ppp[,-ncol(stored_centers_all_ppp)])
  } 
  if (run_type == "osun" | run_type == "known" | run_type == "blm") {
    center_comparison <- compare_centers(truth_ppp = equidistant_truth_list_ppp, centers_df = stored_centers_all_ppp)
  }
  
  # Matrix #1 interpretation: y-axis is further from true P% of the time than x-axis. Note: If values of matrix pairs don't sum up to 1, then, for the remaining difference, distances are equal. The smaller pair value is still better.
  # Matrix #2 interpretation: When y-axis is further from true than x-axis, it is, on average, P% further
  
  distances_df <- center_comparison$distances
  angles_df <- center_comparison$angles
  
  if (!file.exists(paste0(sub(paste0(basename(variables_path), "/"), "CC_Distance_PCT.csv", variables_path))) | !file.exists(paste0(paste0(sub(paste0(basename(variables_path), "/"), "CC_Magnitude_PCT.csv", variables_path))))) {
    write.csv(center_comparison$distance_pct, paste0(sub(paste0(basename(variables_path), "/"), "CC_Distance_PCT.csv", variables_path)))
    write.csv(center_comparison$magnitude_pct, paste0(sub(paste0(basename(variables_path), "/"), "CC_Magnitude_PCT.csv", variables_path)))
  } else {
    print("File already exists")
  }
}


# -------- INITIALIZE POINT PATTERN CLASSIFICATION ------- #

# Extract chosen center coordinates and weights
if (center_selection == "median") {
  chosen_center <- data.frame(x = stored_centers_all_ppp$Mdx, y = stored_centers_all_ppp$Mdy)
  pattern_dataset <- list()
  for (ppp in 1:length(ppp_list_equidistant)) {
    pattern_dataset[[ppp]] <- coords(ppp_list_equidistant[[ppp]])
  }
} else if (center_selection == "mean") {
  chosen_center <- data.frame(x = stored_centers_all_ppp$Mx, y = stored_centers_all_ppp$My)
  pattern_dataset <- list()
  for (ppp in 1:length(ppp_list_equidistant)) {
    pattern_dataset[[ppp]] <- coords(ppp_list_equidistant[[ppp]])
  }
} else if (center_selection == "weighted_mean") {
  chosen_center <- data.frame(x = stored_centers_all_ppp$Mwx, y = stored_centers_all_ppp$Mwy)
  pattern_dataset <- stored_full_dataset_all_ppp
} else if (grepl("constrained_weighted_mean", center_selection)) {
  extract_cols <- stored_centers_all_ppp[grepl("C", colnames(stored_centers_all_ppp))]
  extract_cols <- extract_cols[grepl(gsub("\\D", "", center_selection), colnames(extract_cols))]
  chosen_center <- data.frame(x = extract_cols[,1], y = extract_cols[,2])
  pattern_dataset <- constrained_reweights_list[[as.numeric(gsub("\\D", "", center_selection))]]
}


# ---- STANDARD DISTANCE CALCULATION (COMPACTNESS) ----- #

stdist <- list()
for (i in adjustments) {
  if (center_selection == "median" | center_selection == "mean") {
    stdist[[i]] <- standard_distance(pattern_dataset = pattern_dataset,
                                coordinates_colnames = c("x", "y"),
                                chosen_center = chosen_center,
                                center_colnames = c("x", "y"),
                                adjustments = adjustments[i])
  } else {
    stdist[[i]] <- standard_distance(pattern_dataset = pattern_dataset,
                                coordinates_colnames = c("x", "y"),
                                chosen_center = chosen_center,
                                center_colnames = c("x", "y"),
                                weights_colname = "rect",
                                adjustments = adjustments[i])
  }
  print(paste0(i, " / ", length(adjustments)))
}


# ------------ STANDARD DEVIATIONAL ELLIPSE ------------ #

stored_ellipses <- list()
stored_major <- list()
stored_minor <- list()
stored_theta <- list()
for (i in adjustments) {
  if (center_selection == "median" | center_selection == "mean") {
    sde <- standard_dev_ellipse(pattern_dataset = pattern_dataset,
                                coordinates_colnames = c("x", "y"),
                                chosen_center = chosen_center,
                                center_colnames = c("x", "y"),
                                adjustments = adjustments[i])
  } else {
    sde <- standard_dev_ellipse(pattern_dataset = pattern_dataset,
                                coordinates_colnames = c("x", "y"),
                                chosen_center = chosen_center,
                                center_colnames = c("x", "y"),
                                weights_colname = "rect",
                                adjustments = adjustments[i])
  }
  stored_ellipses[[i]] <- sde$ellipses
  stored_major[[i]] <- sde$major
  stored_minor[[i]] <- sde$minor
  stored_theta[[i]] <- sde$theta
  print(paste0(i, " / ", length(adjustments)))
}


# ------------- WITHIN-ELLIPSE PERCENTAGE -------------- #

if (any(grepl("true", colnames(dataset_sorted)) & grepl("lat", colnames(dataset_sorted)))) {
  dataset_subset <- unique(data.frame(samples = dataset_sorted$sampleID, y = dataset_sorted$lat_true, x = dataset_sorted$long_true))
} else if (any(colnames(dataset_sorted) == "x")) {
  dataset_subset <- unique(data.frame(samples = dataset_sorted$sampleID, y = dataset_sorted$y, x = dataset_sorted$x))
}

if (all(is.na(restriction))) {
  # Rapid Version
  truth_in_ellipse_list <- list()
  percent_sdev_list <- list()
  for (i in adjustments) {
    truth_in_ellipse <- c()
    for (index in 1:length(equidistant_truth_list_ppp)) {
      dat <- data.frame(in.ell = as.logical(sp::point.in.polygon(equidistant_truth_list_ppp[[index]]$x, equidistant_truth_list_ppp[[index]]$y, stored_ellipses[[i]][[index]][,1], stored_ellipses[[i]][[index]][,2])))
      truth_in_ellipse <- c(truth_in_ellipse, dat$in.ell)
    }
    percent_sdev_list[[i]] <- sum(truth_in_ellipse) / length(truth_in_ellipse)
  }
} else {
  # Slow Version
  in_vs_out_list <- list()
  percent_sdev_list <- list()
  for (i in adjustments) {
    in_vs_out <- c()
    for (samples in 1:length(equidistant_truth_list_ppp)) {
      if (all(!is.na(equidistant_truth_list_ppp[[samples]]))) {
        to_poly <- stored_ellipses[[i]][[samples]] %>%
          st_as_sf(coords = c("x", "y"), crs = equidistant) %>%
          summarise(geometry = st_combine(geometry)) %>%
          st_cast("POLYGON") 
        io <- st_intersects(to_poly, st_as_sf(coords(equidistant_truth_list_ppp[[samples]]), coords = c("x", "y"), crs = equidistant))[[1]]
        if (length(io) == 1) {
          in_vs_out <- c(in_vs_out, io)
        } else {
          in_vs_out <- c(in_vs_out, 0)
        }
      } else {
        in_vs_out <- c(in_vs_out, NA)
      }
      percent_sdev <- sum(in_vs_out) / length(in_vs_out)
      print(paste0(samples, " / ", length(equidistant_truth_list_ppp), " || ", i, " / ", length(adjustments)))
    }
    in_vs_out_list[[i]] <- in_vs_out
    percent_sdev_list[[i]] <- percent_sdev
  }
  # RESTRICTED IN_VS_OUT
  restricted_in_ellipse <- list()
  for (i in adjustments) {
    restricted_in_ellipse[[i]] <- sum(in_vs_out_list[[i]][dataset_subset$y >= restriction[1] & dataset_subset$y <= restriction[2]]) / length(in_vs_out_list[[i]][dataset_subset$y >= restriction[1] & dataset_subset$y <= restriction[2]])
  }
}


# ----------------- ADDITIONAL OUTPUTS ----------------- #

# Distance and Angle Errors for tested center types
for (col in 1:ncol(distances_df)) {
  for (i in 1:length(distances_df[,col])) {
    distances_df[,col][i] <- distances_df[,col][i] / 1000
  }
}
distances_df$Samples <- unique(dataset_subset$samples)
distances_df <- distances_df[,c(ncol(distances_df), 1:(ncol(distances_df) - 1))]
distances_df_unrestricted <- distances_df
if (all(!is.na(restriction))) {
  distances_df_restricted <- distances_df[dataset_subset$y >= restriction[1] & dataset_subset$y <= restriction[2],]
}

angles_df$Samples <- unique(dataset_subset$samples)
angles_df <- angles_df[,c(ncol(angles_df), 1:(ncol(angles_df) - 1))]
angles_df_unrestricted <- angles_df
if (all(!is.na(restriction))) {
  angles_df_restricted <- angles_df[dataset_subset$y >= restriction[1] & dataset_subset$y <= restriction[2],]
}

# Standard Deviational Errors for 1 - N standard deviations
A <- Map("*", stored_major, stored_minor)
B <- Map("*", A, pi)
C <- Map("/", B, 10^6)

populated_errors_restricted <- data.frame(Samples = NA, Standard_Deviations = NA, Major_Error = NA, Minor_Error = NA, Area_Error = NA)
populated_errors_unrestricted <- data.frame(Samples = NA, Standard_Deviations = NA, Major_Error = NA, Minor_Error = NA, Area_Error = NA)
for (i in adjustments) {
  if (all(!is.na(restriction))) {
    sam <- dataset_subset$samples[dataset_subset$y >= restriction[1] & dataset_subset$y <= restriction[2]]
    smaj <- stored_major[[i]][dataset_subset$y >= restriction[1] & dataset_subset$y <= restriction[2]]
    smin <- stored_minor[[i]][dataset_subset$y >= restriction[1] & dataset_subset$y <= restriction[2]]
    sare <- C[[i]][dataset_subset$y >= restriction[1] & dataset_subset$y <= restriction[2]]
    populated_errors_restricted <- rbind(populated_errors_restricted, data.frame(Samples = sam, Standard_Deviations = i, Major_Error = smaj / 1000, Minor_Error = smin / 1000, Area_Error = sare))
  }
  populated_errors_unrestricted <- rbind(populated_errors_unrestricted, data.frame(Samples = dataset_subset$samples, Standard_Deviations = i, Major_Error = stored_major[[i]] / 1000, Minor_Error = stored_minor[[i]] / 1000, Area_Error = C[[i]]))
}
if (all(!is.na(restriction))) {
  populated_errors_restricted <- populated_errors_restricted[2:nrow(populated_errors_restricted),]
}
populated_errors_unrestricted <- populated_errors_unrestricted[2:nrow(populated_errors_unrestricted),]

write.csv(distances_df_unrestricted, paste0(sub(paste0(basename(variables_path), "/"), "DISTANCE_ERRORS.csv", variables_path)))
write.csv(angles_df_unrestricted, paste0(sub(paste0(basename(variables_path), "/"), "ANGLE_OFFSET.csv", variables_path)))
write.csv(populated_errors_unrestricted, paste0(sub(paste0(basename(variables_path), "/"), "ELLIPSE_ERRORS.csv", variables_path)))

if (all(!is.na(restriction))) {
  write.csv(distances_df_restricted, paste0(sub(paste0(basename(variables_path), "/"), "RESTRICTED_DISTANCE_ERRORS.csv", variables_path)))
  write.csv(angles_df_restricted, paste0(sub(paste0(basename(variables_path), "/"), "RESTRICTED_ANGLE_OFFSET.csv", variables_path)))
  write.csv(populated_errors_restricted, paste0(sub(paste0(basename(variables_path), "/"), "RESTRICTED_ELLIPSE_ERRORS.csv", variables_path)))
}

# ---------- INDIVIDUAL SAMPLE VISUALIZATION ----------- #

if (all(!is.na(restriction))) {
  dataset_which <- which(dataset_subset$y >= 43 & dataset_subset$y <= 47)
} else {
  dataset_which <- 1:length(unique(dataset_subset$samples))
}

if (!file.exists(paste0(visuals_path, center_selection, "_Distances_Plots.pdf"))) {
  pdf(paste0(visuals_path, center_selection, "_Distances_Plots.pdf"))
  n <- 0
  for (index in dataset_which) {
    n <- n + 1
    # plot(study_area_equalarea, xlim = c(-476455.0, 339383.8), ylim = c(-601974.6, 1007669.4))
    plot(study_area_equidistant)
    points(ppp_list_equidistant[[index]]$x, ppp_list_equidistant[[index]]$y, pch = 16, col = "lightgrey", cex = 0.5)
    for (i in adjustments) {
      lines(stored_ellipses[[i]][[index]], col = "red", lwd = 1)
    }
    points(chosen_center[index,]$x, chosen_center[index,]$y, bg = "purple", pch = 21, lwd = 1, cex = 1.2)
    if (truth.exists) {
      if (any(!is.na(equidistant_truth_list_ppp[[index]]))) {
        points(equidistant_truth_list_ppp[[index]]$x, equidistant_truth_list_ppp[[index]]$y, bg = "green", pch = 21, lwd = 1, cex = 1.2)
      }
    }
    if (all(!is.na(restriction))) {
      title(main = paste0(distances_df_restricted$Samples[n], " | Error: ", distances_df_restricted$Distance_WMN[n], "km"), cex.main = 0.7)
    } else {
      title(main = paste0(distances_df_unrestricted$Samples[n], " | Error: ", distances_df_unrestricted$Distance_WMN[n], "km"), cex.main = 0.7)
    }
  }
  dev.off()
} else {
  print("File Already Exists")
}
