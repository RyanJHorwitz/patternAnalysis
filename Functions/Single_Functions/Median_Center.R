
# Median_Center.R

# Ported from Python (PySAL) https://pysal.org/notebooks/explore/pointpats/centrography.html
# pts: A spatstat ppp object
# crit: The convergence criterion

median_center <- function(pts, crit = 0.0001) {
  x <- pts[,1]
  y <- pts[,2]
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
  return (output)
}
