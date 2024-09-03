
# Standard_Distance.R

# Source: https://pro.arcgis.com/en/pro-app/latest/tool-reference/spatial-statistics/h-how-standard-distance-spatial-statistic-works.htm

# pattern_dataset: A list of spatstat ppp objects
# coordinates_colnames: Vector of column names (as string) of the X and Y coordinates, respectively
# chosen_center: A spatstat ppp object containing the XY coordinates representing the center of the point pattern. For instance, mean center, median center, etc.
# weights_colname: The column name of the weights column (if applicable)
# adjustments: number of (weighted) standard distances/deviations from the (weighted) mean

standard_distance <- function(pattern_dataset, coordinates_colnames, chosen_center, center_colnames, weights_colname = NA, adjustments = 1) {

  # Weighted Standard Distance
  if (!is.na(weights_colname)) {
    SDw <- c()
    for (index in 1:length(pattern_dataset)) {
      top_x <- c()
      top_y <- c()
      bottom <- c()
      for (point in 1:nrow(pattern_dataset[[index]])) {
        wi <- pattern_dataset[[index]][,colnames(pattern_dataset[[index]]) == weights_colname][point] # weight of point
        xi <-   pattern_dataset[[index]][,colnames(pattern_dataset[[index]]) == coordinates_colnames[1]][point] # x-coordinate of point
        xw <-  chosen_center[index,][,colnames(chosen_center[index,]) == center_colnames[1]] # x-coordinate of center
        yi <-   pattern_dataset[[index]][,colnames(pattern_dataset[[index]]) == coordinates_colnames[2]][point] # y-coordinate of point
        yw <- chosen_center[index,][,colnames(chosen_center[index,]) == center_colnames[2]] # y-coordinate of center
        top_x <- c(top_x, wi*((xi - xw)^2))
        top_y <- c(top_y, wi*((yi - yw)^2))
        bottom <- c(bottom, wi)
      }
      SDw <- c(SDw, sqrt((sum(top_x) / sum(bottom)) + (sum(top_y) / sum(bottom))))
    }
    SDw <- SDw * adjustments
    output <- SDw
  } 
  # Unweighted Standard Distance
  else {
    SD <- c()
    for (index in 1:length(pattern_dataset)) {
      top_x <- c()
      top_y <- c()
      bottom <- nrow(pattern_dataset[[index]])
      for (point in 1:nrow(pattern_dataset[[index]])) {
        xi <-   pattern_dataset[[index]][,colnames(pattern_dataset[[index]]) == coordinates_colnames[1]][point] # x-coordinate of point
        xw <- chosen_center$x[index] # x-coordinate of center
        yi <-   pattern_dataset[[index]][,colnames(pattern_dataset[[index]]) == coordinates_colnames[2]][point] # y-coordinate of point
        yw <- chosen_center$y[index] # y-coordinate of center
        top_x <- c(top_x, ((xi - xw)^2))
        top_y <- c(top_y, ((yi - yw)^2))
      }
      SD <- c(SD, sqrt((sum(top_x) / bottom) + (sum(top_y) / bottom)))
    }
    SD <- SD * adjustments
    output <- SD
  }
  return(output)
}
