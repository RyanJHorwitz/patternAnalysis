
# Mean_Center.R

# Take the average of all x/y coordinates for each elephant point pattern

# ppp_object: A spatstat ppp object

mean_center <- function(ppp_object) {

  mean_x <- mean(coords(ppp_object)[,1]) # Mean x-coordinate
  mean_y <- mean(coords(ppp_object)[,2]) # Mean y-coordinate

  output <- c(x = mean_x, y = mean_y)
  
  return(output)
  
}