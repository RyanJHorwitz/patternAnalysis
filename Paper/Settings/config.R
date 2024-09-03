
# FILE: config.R

# ----------------------- OPTIONS ----------------------- #

species <- ''
run_type <- ''
sim_version <- NA
drop_version <- NA
data_selection <- ''
study_area <- ''
equalarea <- ''
equidistant <- ''
adjustments <- c()
constraint_paths <- c()
center_selection <- ''
restriction <- NA
neighborhood_size <- 0.5

# species: Options ... "African_Forest_Elephant", "African_Savannah_Elephant", "Douglas_Fir", "Big_Leaf_Maple"
# run_type: Options ... (elephant: "Simulated", "Real", "Drop"), (douglas fir: "Known", "OSUN"), (Big Leaf Maple: "BLM")
# sim_version: Options ... Only run_type "Simulated": String specifying version of simulated data (typically by date)
# drop_version: Options ... Only run_type "Drop": String specifying type of drop data (e.g. "Jackknife_Predictions)
# data_selection: Options ... Always set to "All" unless run_type "Real", in which case can set to all or "Subset.csv" file
# study_area: Options ... Path to polygon shapefile  containing the study area to contextualize the point data. Used in output maps
# equalarea: Options ... Map projection that minimizes area distortion of the dataset
# equidistant: Options ... Map projection that minimizes distance distortion of the dataset
# adjustments: Options ... Number of standard deviations to calculate. Can be integer or vector of integers
# constraint_paths: Options ... Path to polygon shapefile or vector of paths to shapefiles representing the known habitat of the species of interest. Often with buffer. Used in constrained weighted mean center calculation
# center_selection: Options ... FALSE or "Mean", "Median", "Weighted_Mean, "Constrained_Weighted_Mean_[n]" (where [n] corresponds to integer ID of constraining buffer)
# restriction: Options ... NA or vector of 2 latitudes between which the data sits. Removes any points outside of those latitudes if set.
# neighborhood_size: Options ... Decimal between 0 and 1 representing proportion of samples to use (percent) or integer if 1 or greater representing number of samples (count)


# --------------------------------------------------------------------- #
# --------------------------------------------------------------------- #
# --- IGNORE --- IGNORE --- IGNORE --- IGNORE --- IGNORE --- IGNORE --- #
# --- IGNORE --- IGNORE --- IGNORE --- IGNORE --- IGNORE --- IGNORE --- #
# --- IGNORE --- IGNORE --- IGNORE --- IGNORE --- IGNORE --- IGNORE --- #
# --- IGNORE --- IGNORE --- IGNORE --- IGNORE --- IGNORE --- IGNORE --- #
# --- IGNORE --- IGNORE --- IGNORE --- IGNORE --- IGNORE --- IGNORE --- #
# --- IGNORE --- IGNORE --- IGNORE --- IGNORE --- IGNORE --- IGNORE --- #
# --------------------------------------------------------------------- #
# --------------------------------------------------------------------- #

# ----------------------------------------------------- #
# ------------------ POST-PROCESSING ------------------ #
# ------------------ POST-PROCESSING ------------------ #
# ------------------ POST-PROCESSING ------------------ #
# ------------------ POST-PROCESSING ------------------ #
# ----------------------------------------------------- #


# ------------- SETTINGS INPUT PROCESSING ------------- #

# ---- INPUT MODIFICATIONS

# species
if (species == "African_Forest_Elephant") {
  constraint_paths <- str_replace(constraint_paths, species, "forest")
  species <- "forest"
}
if (species == "African_Savannah_Elephant") {
  constraint_paths <- str_replace(constraint_paths, species, "savannah")
  species <- "savannah"
}
if (species == "Douglas_Fir") {
  constraint_paths <- str_replace(constraint_paths, species, "douglas_fir")
  species <- "douglas_fir"
}
if (species == "Big_Leaf_Maple") {
  constraint_paths <- str_replace(constraint_paths, species, "big_leaf_maple")
  species <- "big_leaf_maple"
}

# constraint_paths
constraint_paths <- paste0('./Settings/dependencies/', constraint_paths)

# run_type
run_type <- tolower(run_type)
if (any(run_type == "simulated" | run_type == "real" | run_type == "drop")) {
  split.dat <- T
  iterations <- 9
}
if (any(run_type == "known" | run_type == "osun" | run_type == "blm")) {
  split.dat <- F
}

# data_selection
data_selection <- tolower(data_selection)
if (data_selection == "subset.csv") {
  data_selection <- "./Settings/dependencies/subset.csv"
  data_selection <- c(read.csv(data_selection, header = F, fileEncoding="UTF-8-BOM"))[[1]]
}

# study_area
study_area <- paste0('./Settings/dependencies/', study_area)
study_area <- st_geometry(read_sf(dsn = study_area)) # Read in data
study_area_equalarea <- st_transform(study_area, crs = equalarea) # Project Coordinates
study_area_equidistant <- st_transform(study_area, crs = equidistant) # Project Coordinates

# center_selection
if (center_selection != FALSE) {
  center_selection <- tolower(center_selection)
}

# ---- ERROR CALLS

if (any(run_type == "simulated" | run_type == "real" | run_type == "drop") & 
    (!any(species == "forest" | species == "savannah"))) {
  stop("run_type not compatible with species")
}

if (any(run_type == "known" | run_type == "osun") & 
    (!any(species == "douglas_fir"))) {
  stop("run_type not compatible with species")
}

if (run_type == "simulated" & sim_version == "") {
  stop("run_type simulated requires sim_version variable to be specified")
}

if (run_type == "drop" & drop_version == "") {
  stop("run_type drop requires drop_version variable to be specified")
}

if ((run_type != "real") & (!any(data_selection == "all"))) {
  stop("data_selection not compatible with run_type")
}

if (any(typeof(adjustments) != "integer") | any(is.na(adjustments))) {
  stop("adjustments must be integer")
}


# ------- PREPARING ANALYSIS INPUTS AND OUTPUTS ------- #

# ---- run_type set to REAL
if (run_type == "real") {
  input_path <- "./Datasets/elephant/real"
  if (!dir.exists(input_path)) {
    stop("Path to dataset does not exist")
  }
  input_list <- list.dirs(input_path, recursive = F)
  if (data_selection == "all") {
    variables_path <- paste0("./Outputs/", run_type, "/", species, "/full_dataset/variables/")
    visuals_path <- paste0("./Outputs/", run_type, "/", species, "/full_dataset/visuals/")
  } else {
    input_list <- list.dirs(input_path, recursive = F)
    input_list <- input_list[basename(input_list) %in% data_selection]
    variables_path <- paste0("./Outputs/", run_type, "/", species, "/select_seizure_runs/", paste(data_selection, collapse = "-"), "/variables/")
    visuals_path <- paste0("./Outputs/", run_type, "/", species, "/select_seizure_runs/", paste(data_selection, collapse = "-"), "/visuals/")
  }
  truth.exists <- FALSE
}

# ---- run_type set to SIMULATED
if (run_type == "simulated") {
  input_path <- paste0("./Datasets/elephant/", run_type, "/", sim_version, "/", species, "/")
  if (!dir.exists(input_path)) {
    stop("Path to dataset does not exist")
  }
  input_list <- list.dirs(input_path, recursive = F)
  input_list <- input_list[which(!grepl("true", input_list))]
  variables_path <- paste0("./Outputs/", run_type, "/", sim_version, "/", species, "/variables/")
  visuals_path <- paste0("./Outputs/", run_type, "/", sim_version, "/", species, "/visuals/")
  truth.exists <- TRUE
  true_path <- list.dirs(input_path, recursive = F)
  true_path <- true_path[which(grepl("true", true_path))]
}

# ---- run_type set to DROP
if (run_type == "drop") {
  input_path <- paste0("./Datasets/elephant/", run_type, "/", drop_version, "/predictions/", species, "/")
  if (!dir.exists(input_path)) {
    stop("Path to dataset does not exist")
  }
  variables_path <- paste0("./Outputs/", run_type, "/", drop_version, "/", species, "/variables/")
  visuals_path <- paste0("./Outputs/", run_type, "/", drop_version, "/", species,  "/visuals/")
  truth.exists <- TRUE
  true_path <- list.files(paste0("./Datasets/elephant/", run_type, "/"), full.names = T)
  true_path <- true_path[grepl("true", true_path)]
}

# ---- run_type set to KNOWN or OSUN
if (run_type == "known" | run_type == "osun") {
  input_path <- paste0("./Datasets/douglas_fir/", run_type, "/")
  if (!dir.exists(input_path)) {
    stop("Path to dataset does not exist")
  }
  variables_path <- paste0("./Outputs/", run_type, "/variables/")
  visuals_path <- paste0("./Outputs/", run_type, "/visuals/")
  truth.exists <- TRUE
}

# ---- run_type set to BLM
if (run_type == "blm") {
  input_path <- paste0("./Datasets/big_leaf_maple/", run_type, "/")
  if (!dir.exists(input_path)) {
    stop("Path to dataset does not exist")
  }
  variables_path <- paste0("./Outputs/", run_type, "/variables/")
  visuals_path <- paste0("./Outputs/", run_type, "/visuals/")
  truth.exists <- TRUE
}

# ---- CREATE OUTPUT DIRS IF DON'T EXIST

if (!dir.exists("../Outputs")) {
  dir.create("../Outputs")
}
if (!dir.exists(variables_path)) {
  dir.create(variables_path, recursive = T)
}
if (!dir.exists(visuals_path)) {
  dir.create(visuals_path, recursive = T)
}

# --------- READ-IN SCAT PREDICTION DATAFRAMES --------- #

# --- run type is REAL or SIMULATED
if (run_type == "real" | run_type == "simulated") {
  data_fullpaths <- c()
  stored_group_lengths <- c()
  i <- 0
  for (folder in input_list) {
    i <- i + 1
    estimates_folder <- list.dirs(folder, recursive = F)
    iterations_folder <- list.dirs(estimates_folder[grepl(paste0("n", species), estimates_folder)], recursive = F)
    outputs_folder <- list.dirs(iterations_folder[nchar(basename(iterations_folder)) == 1], recursive = F)
    data_files <- unique(list.files(outputs_folder))
    data_names <- data_files[!grepl("Output", data_files)]
    data_fullpaths <- c(data_fullpaths, paste0(rep(paste0(outputs_folder, "/"), length(data_names)), sort(rep(data_names, length(outputs_folder)))))
    stored_group_lengths <- c(stored_group_lengths, length(data_names))
  }
  extract <- seq(from = 1, to = length(data_fullpaths), by = iterations)
  all_data_names <- basename(data_fullpaths)[extract]
  if (!file.exists(paste0(variables_path, "predictions.RDS"))) {
    i <- 0
    n <- 0
    m <- 0
    sample_group <- data.frame()
    predictions <- list()
    for (tables in data_fullpaths) {
      i <- i + 1
      m <- m + 1
      if ((i > iterations)) {
        n <- n + 1
        predictions[[n]] <- sample_group
        sample_group <- data.frame()
        i <- 1
      }
      data_table <- read.table(tables)
      data_table <- data_table[(nrow(data_table) - 100):(nrow(data_table) - 1), 1:2]
      sample_group <- data.frame(rbind(sample_group, data_table))
      rownames(sample_group) <- NULL
      if (m == length(data_fullpaths)) {
        n <- n + 1
        predictions[[n]] <- sample_group
      }
      print(paste0(m, " / ", length(data_fullpaths)))
    }
    saveRDS(predictions, paste0(variables_path, "predictions.RDS"))
  } else {
    predictions <- readRDS(paste0(variables_path, "predictions.RDS"))
  }
}

# --- run_type is DROP
if (run_type == "drop") {
  data_fullpaths <- list.files(input_path, recursive = F, full.names = T)
  data_fullpaths <- data_fullpaths[!str_detect(data_fullpaths, ".txt")]
  stored_group_lengths <- length(data_fullpaths)
  all_data_names <- sub("^(([^_]*_){2}[^_]*).*", "\\1", basename(data_fullpaths))
  drop_zones <- sub("^(([^_]*_){1}[^_]*).*", "\\1", all_data_names)
  drop_lengths <- data.frame(data.frame(drop_zones) %>% group_by(drop_zones) %>% summarise(total_count = n()))
  if (!file.exists(paste0(variables_path, "predictions.RDS"))) {
    predictions <- list()
    n <- 0
    for (tables in data_fullpaths) {
      n <- n + 1
      data_table <- read.table(tables)
      predictions[[n]] <- data.frame(V1 = data_table$V2, V2 = data_table$V1)
      print(paste0(n, " / ", length(data_fullpaths)))
    }
    saveRDS(predictions, paste0(variables_path, "predictions.RDS"))
  } else {
    predictions <- readRDS(paste0(variables_path, "predictions.RDS"))
  }
}

# --- run_type is KNOWN or OSUN
if (run_type == "known" | run_type == "osun") {
  data_fullpaths <- list.files(input_path, recursive = F, full.names = T)
  dataset <- read.csv(data_fullpaths)
  dataset_sorted <- dataset[order(dataset$sampleID),]
  all_data_names <- unique(dataset_sorted$sampleID)
  stored_sample_lengths <- c(table(dataset_sorted$sampleID))
  if (!file.exists(paste0(variables_path, "predictions.RDS"))) {
    tables_list <- split(dataset_sorted , f = dataset_sorted$sampleID)
    predictions <- list()
    n <- 0
    for (tables in tables_list) {
      n <- n + 1
      predictions[[n]] <- data.frame(V1 = tables$lat_est, V2 = tables$long_est)
      print(paste0(n, " / ", length(tables_list)))
    }
    saveRDS(predictions, paste0(variables_path, "predictions.RDS"))
  } else {
    predictions <- readRDS(paste0(variables_path, "predictions.RDS"))
  }
}

# --- run_type is BLM
if (run_type == "blm") {
  data_fullpaths <- list.files(input_path, recursive = F, full.names = T)
  dataset <- read.csv(data_fullpaths[1])
  dataset_sorted <- dataset[order(dataset$sampleID),]
  all_data_names <- unique(dataset_sorted$sampleID)
  stored_sample_lengths <- c(table(dataset_sorted$sampleID))
  if (!file.exists(paste0(variables_path, "predictions.RDS"))) {
    tables_list <- split(dataset_sorted , f = dataset_sorted$sampleID)
    predictions <- list()
    n <- 0
    for (tables in tables_list) {
      n <- n + 1
      predictions[[n]] <- data.frame(V1 = tables$y, V2 = tables$x)
      print(paste0(n, " / ", length(tables_list)))
    }
    saveRDS(predictions, paste0(variables_path, "predictions.RDS"))
  } else {
    predictions <- readRDS(paste0(variables_path, "predictions.RDS"))
  }
}

# --------- CONVERSION TO PROJECTED SF OBJECTS --------- #

if (!file.exists(paste0(variables_path, "locations_list_equidistant.RDS")) | !file.exists(paste0(variables_path, "locations_list_equalarea.RDS"))) {
  locations_list_equidistant <- list()
  locations_list_equalarea <- list()
  for (index in 1:length(predictions)) {
    spatial_locations <- st_as_sf(predictions[[index]], coords = c("V2","V1")) # Convert dataframe to simple feature
    spatial_locations <- st_set_crs(spatial_locations, 4326) # Set geographic coordinate system of the point data
    locations_projected_equidistant <- st_transform(spatial_locations, crs = equidistant) # Project the coordinates of the point data
    locations_projected_equalarea <- st_transform(spatial_locations, crs = equalarea) # Project the coordinates of the point data
    locations_list_equidistant[[index]] <- locations_projected_equidistant # Store simple feature for sample in list with all other samples
    locations_list_equalarea[[index]] <- locations_projected_equalarea # Store simple feature for sample in list with all other samples
    print(paste0(index, " / ", length(predictions)))
  }
  saveRDS(locations_list_equidistant, paste0(variables_path, "locations_list_equidistant.RDS"))
  saveRDS(locations_list_equalarea, paste0(variables_path, "locations_list_equalarea.RDS"))
} else {
  locations_list_equidistant <- readRDS(paste0(variables_path, "locations_list_equidistant.RDS"))
  locations_list_equalarea <- readRDS(paste0(variables_path, "locations_list_equalarea.RDS"))
  print("SF Objects Successfully Loaded")
}

# -------------- CONVERSION TO PPP OBJECTS ------------- #

if (!file.exists(paste0(variables_path, "ppp_list_equidistant.RDS")) | !file.exists(paste0(variables_path, "ppp_list_equalarea.RDS"))) {
  ppp_list_equidistant <- list()
  ppp_list_equalarea <- list()
  for (i in 1:length(locations_list_equidistant)) {
    ppp_list_equidistant[[i]] <- as.ppp(st_geometry(locations_list_equidistant[[i]]), as.owin(st_bbox(locations_list_equidistant[[i]])))
    ppp_list_equalarea[[i]] <- as.ppp(st_geometry(locations_list_equalarea[[i]]), as.owin(st_bbox(locations_list_equalarea[[i]])))
    print(paste0(i, " / ", length(locations_list_equidistant)))
  }
  saveRDS(ppp_list_equidistant, paste0(variables_path, "ppp_list_equidistant.RDS"))
  saveRDS(ppp_list_equalarea, paste0(variables_path, "ppp_list_equalarea.RDS"))
} else {
  ppp_list_equidistant <- readRDS(paste0(variables_path, "ppp_list_equidistant.RDS"))
  ppp_list_equalarea <- readRDS(paste0(variables_path, "ppp_list_equalarea.RDS"))
  print("PPP Data Loaded from File")
}

# ------ IF TRUTH EXISTS, LOAD IN TRUE COORDINATES ----- #

if (truth.exists) {
  
  # If run type is "SIMULATED"
  if (run_type == "simulated") {
    truth <- list.files(true_path, full.names = T)
    truth_df <- data.frame()
    for (table in truth) {
      truth_df <- data.frame(rbind(truth_df, read.table(table, header = T)))
    }
    colnames(truth_df) <- c("Elephant", "zone_y", "zone_x", "individual_y", "individual_x")
    truth_df$zone_y <- as.numeric(gsub(",","", truth_df$zone_y))
    truth_df$individual_y <- as.numeric(gsub(",","", truth_df$individual_y))
    truth_df <- data.frame(cbind(data.frame(Group = rep(basename(input_list), stored_group_lengths)), truth_df))
    truth_sf <- st_as_sf(truth_df, coords = c("zone_x", "zone_y"), crs = 4326)
    equalarea_truth_transformed <- st_transform(truth_sf, crs = equalarea)
    equidistant_truth_transformed <- st_transform(truth_sf, crs = equidistant)
    equalarea_truth_list <- list()
    equidistant_truth_list <- list()
    for (points in 1:nrow(equalarea_truth_transformed)) {
      equalarea_truth_list[[points]] <- equalarea_truth_transformed[points,]
      equidistant_truth_list[[points]] <- equidistant_truth_transformed[points,]
    }
    equalarea_truth_list_ppp <- list()
    equidistant_truth_list_ppp <- list()
    for (index in 1:length(equalarea_truth_list)) {
      equalarea_truth_list_ppp[[index]] <- as.ppp(st_geometry(equalarea_truth_list[[index]]), as.owin(st_bbox(equalarea_truth_list[[index]])))
      equidistant_truth_list_ppp[[index]] <- as.ppp(st_geometry(equidistant_truth_list[[index]]), as.owin(st_bbox(equidistant_truth_list[[index]])))
      print(paste0(index, " / ", length(equalarea_truth_list)))
    }
  }
  
  # If run type is "DROP"
  if (run_type == "drop") {
    truth <- true_path
    truth_df <- read.csv(truth, header = T, row.names = NULL)
    truth_df <- truth_df[colnames(truth_df) %in% c("pred_file", "Input.Zone", "Match.ID", "Species", "Sample.ID", "Loci", "Region", "Sector", "Sector.Name", "Location", "Latitude", "Longitude", "Country", "Sample.type", "Zone")]
    colnames(truth_df) <- c("pred_file", "input_zone", "match_id", "species", "sample_id", "loci", "region", "sector", "sector_name", "location", "lat", "long", "country", "sample_type", "zone")
    truth_df <- truth_df[truth_df$species == toupper(substr(species, 1, 1)),]
    truth_df <- truth_df[-which(is.na(truth_df$long) | is.na(truth_df$lat)),]
    
    subset_truth_names <- sub("_predictions.tsv", "", truth_df$pred_file)
    ordered_truth_names <- subset_truth_names[order(match(subset_truth_names, all_data_names))]
    ordered_truth_df <- truth_df[order(match(subset_truth_names, all_data_names)),]
    truth_df_list <- split(ordered_truth_df, f = factor(paste0(ordered_truth_names, "_predictions.tsv"), levels = unique(paste0(ordered_truth_names, "_predictions.tsv"))))
    truth_sf <- c()
    i <- 0
    for (dfs in truth_df_list) {
      i <- i + 1
      dfs <- dfs[!duplicated(cbind(dfs$long, dfs$lat)),]
      center_sf <- st_as_sf(dfs, coords = c("long", "lat"), crs = 4326)
      if (nrow(center_sf) > 1) {
        single_center <- median_center(st_coordinates(center_sf))
        dfs$long <- single_center[1]
        dfs$lat <- single_center[2]
        dfs <- dfs[!duplicated(cbind(dfs$long, dfs$lat)),]
        center_sf <- st_as_sf(dfs, coords = c("long", "lat"), crs = 4326)
      }
      truth_sf <- rbind(truth_sf, center_sf)
    }
    equalarea_truth_transformed <- st_transform(truth_sf, crs = equalarea)
    equidistant_truth_transformed <- st_transform(truth_sf, crs = equidistant)
    equalarea_truth_list <- list()
    equidistant_truth_list <- list()
    for (points in 1:nrow(equalarea_truth_transformed)) {
      equalarea_truth_list[[points]] <- equalarea_truth_transformed[points,]
      equidistant_truth_list[[points]] <- equidistant_truth_transformed[points,]
    }
    equalarea_truth_list_ppp <- list()
    equidistant_truth_list_ppp <- list()
    for (index in 1:length(equalarea_truth_list)) {
      equalarea_truth_list_ppp[[index]] <- as.ppp(st_geometry(equalarea_truth_list[[index]]), as.owin(st_bbox(equalarea_truth_list[[index]])))
      equidistant_truth_list_ppp[[index]] <- as.ppp(st_geometry(equidistant_truth_list[[index]]), as.owin(st_bbox(equidistant_truth_list[[index]])))
      print(paste0(index, " / ", length(equalarea_truth_list)))
    }
    equalarea_truth_transformed_df <- data.frame(equalarea_truth_transformed)[,-ncol(equalarea_truth_transformed)]
    equidistant_truth_transformed_df <- data.frame(equidistant_truth_transformed)[,-ncol(equidistant_truth_transformed)]
    equalarea_truth_transformed_df$pred_file <- sub("_[^_]+$", "", equalarea_truth_transformed_df$pred_file)
    equidistant_truth_transformed_df$pred_file <- sub("_[^_]+$", "", equidistant_truth_transformed_df$pred_file)
    identify_no_truth <- left_join(data.frame(samples = all_data_names), equalarea_truth_transformed_df, by = c("samples" = "pred_file"))
    for (index in (which(is.na(identify_no_truth$match_id)) - 1)) {
      equalarea_truth_list <- append(equalarea_truth_list, NA, after = index + 1)
      equalarea_truth_list_ppp <- append(equalarea_truth_list_ppp, NA, after = index + 1)
      equidistant_truth_list <- append(equidistant_truth_list, NA, after = index + 1)
      equidistant_truth_list_ppp <- append(equidistant_truth_list_ppp, NA, after = index + 1)
    }
  }
  
  # If run type is "KNOWN"
  if (run_type == "known") {
    truth <- dataset_sorted
    true_list <- split(dataset_sorted , f = dataset_sorted$sampleID)
    equalarea_truth_list_ppp <- list()
    equidistant_truth_list_ppp <- list()
    n <- 0
    for (tables in true_list) {
      n <- n + 1
      true_points <- data.frame(V1 = tables$lat_true, V2 = tables$long_true)
      truth_sf <- st_as_sf(unique(true_points), coords = c("V2", "V1"), crs = 4326)
      equalarea_truth_list_ppp[[n]] <- as.ppp(st_geometry(st_transform(truth_sf, crs = equalarea)), as.owin(st_bbox(st_transform(truth_sf, crs = equalarea))))
      equidistant_truth_list_ppp[[n]] <- as.ppp(st_geometry(st_transform(truth_sf, crs = equidistant)), as.owin(st_bbox(st_transform(truth_sf, crs = equidistant))))
      print(paste0(n, " / ", length(true_list)))
    }
  }
  
  # If run type is "OSUN"
  if (run_type == "osun") {
    truth <- dataset_sorted
    true_list <- split(dataset_sorted , f = dataset_sorted$sampleID)
    equalarea_truth_list_ppp <- list()
    equidistant_truth_list_ppp <- list()
    n <- 0
    for (tables in true_list) {
      n <- n + 1
      true_points <- data.frame(V1 = tables$lat_true, V2 = tables$long_true)
      if (!any(is.na(true_points))) {
        truth_sf <- st_as_sf(unique(true_points), coords = c("V2", "V1"), crs = 4326)
        equalarea_truth_list_ppp[[n]] <- as.ppp(st_geometry(st_transform(truth_sf, crs = equalarea)), as.owin(st_bbox(st_transform(truth_sf, crs = equalarea))))
        equidistant_truth_list_ppp[[n]] <- as.ppp(st_geometry(st_transform(truth_sf, crs = equidistant)), as.owin(st_bbox(st_transform(truth_sf, crs = equidistant))))
      } else {
        equalarea_truth_list_ppp[[n]] <- NA
        equidistant_truth_list_ppp[[n]] <- NA
      }
      print(paste0(n, " / ", length(true_list)))
    }
  }
  
  # If run type is "BLM"
  if (run_type == "blm") {
    truth <- dataset_sorted
    true_table <- read.table(list.files(input_path, recursive = F, full.names = T)[2], header = T)
    joined_truth <- left_join(truth, true_table, by = "sampleID")
    colnames(joined_truth) <- c("sample_x", "sample_y", "sampleID", "true_x", "true_y")
    true_list <- split(joined_truth , f = joined_truth$sampleID)
    equalarea_truth_list_ppp <- list()
    equidistant_truth_list_ppp <- list()
    n <- 0
    for (tables in true_list) {
      n <- n + 1
      true_points <- data.frame(V1 = tables$true_y, V2 = tables$true_x)
      if (!any(is.na(true_points))) {
        truth_sf <- st_as_sf(unique(true_points), coords = c("V2", "V1"), crs = 4326)
        st_crs(truth_sf) <- 4326
        equalarea_truth_list_ppp[[n]] <- as.ppp(st_geometry(st_transform(truth_sf, crs = equalarea)), as.owin(st_bbox(st_transform(truth_sf, crs = equalarea))))
        equidistant_truth_list_ppp[[n]] <- as.ppp(st_geometry(st_transform(truth_sf, crs = equidistant)), as.owin(st_bbox(st_transform(truth_sf, crs = equidistant))))
      } else {
        equalarea_truth_list_ppp[[n]] <- NA
        equidistant_truth_list_ppp[[n]] <- NA
      }
      print(paste0(n, " / ", length(true_list)))
    }
  }
  
}
