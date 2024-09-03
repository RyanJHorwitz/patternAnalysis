# Point Pattern Analysis

This work proposes methods to best classify and summarize point patterns representing predictions of floral/faunal sample geographic origin.

### Directories

1. Functions: Contains single and batch functions for all measures of central tendency listed below. The batch functions also include a comparison function that suggests the best measure of central tendency for a pattern dataset based on prediction distance to true.

2. Paper: Code used in paper (not yet pubished) to classify geographic origin predictions of douglas fir and big leaf maple samples. The code is sanitized and not all dependencies are included. However, it can serve as a decent example of how to use the functions in the "Functions" folder.

### Point Pattern Central Tendency

The central tendency measure of a point pattern is a single point that attempts to summarize an entire pattern.  Numerous two-dimensional measures of central tendency exist. Among them...

1. The median center minimizes the straight-line distance to all points within the input pattern.  

2. The mean center averages the x and y-coordinates of all points within the input pattern.  

3. The weighted mean center is calculated similarly to the mean center, except each point in the pattern is assigned a weight based on some variable.  In the case of this work, for a single sample point pattern, points are assigned weights based on local density.  The local density of each point, P, is calculated by drawing a minimum bounding rectangle around P and the nearest N points.  The nearest N points + 1 (to account for P) is then divided by the area of the bounding rectangle, resulting in a local density value that serves as the weight of P.

4. The constrained weighted mean is calculated similarly to the weighted mean, except that points within each sample pattern are assigned weights of 0 if they fall outside of a buffered species range of each species of interest.

### Point Pattern Error

The error of a point pattern characterizes the spread of the pattern around the measure of central tendency. The standard deviational ellipse is an ideal measure to classify error, capturing two axes of error and the direction of error.

### Sources

1. Median Center adapted from: https://pysal.org/notebooks/explore/pointpats/centrography.html

2. (Weighted) Mean Center adapted from: https://pro.arcgis.com/en/pro-app/latest/tool-reference/spatial-statistics/mean-center.htm

3. Standard Deviational Ellipse adapted from: https://pro.arcgis.com/en/pro-app/latest/tool-reference/spatial-statistics/h-how-directional-distribution-standard-deviationa.htm
