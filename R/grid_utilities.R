

#######################################
# Grid Build Methods
#######################################

#' @description Build a grid of randomly distribute boxes.
#'     It first creates the stratum, an evenly spaced grid with the defined.
#'     xlimits and ylimits, and then within each stratum grid cell randomly
#'     places boxes. These random boxes are the ones used in the model, and
#'     the stratum grid is discarded after this step.
#'
#' @param stratum_size_x numeric
#' @param stratum_size_y numeric
#' @param boxes_per_stratum numeric
#' @param box_size numeric
#' @param xlimits vector of 2 numerics
#' @param ylimits vector of 2 numerics
#'
#' @return data.frame of box descriptions (centers and sizes)
create_grid_boxes = function(stratum_size_x=0.1,
                       stratum_size_y=0.1,
                       boxes_per_stratum=5,
                       box_size=0.3,
                       xlimits=c(0,1),
                       ylimits=c(0,1)){
  #TODO: add checks for reasonable numbers, densities of boxes.
  # things could get outa hand with a wide limits (ie c(0,100)) and low
  # stratum size (ie. 0.1)

  # Sanity Checks
  if((xlimits[2]<xlimits[1]) |
     (ylimits[2]<ylimits[1]) |
     (length(xlimits) !=2  ) |
     (length(ylimits) !=2  )){stop(paste('xlimits and ylimits must each be of length 2, and have the lower',
                                         'and upper bounds in the 1st and 2nd nposition, respectively.'))}

  # This data.frame describes the bottom left coordinate of the stratum cells.
  stratum_cells = expand.grid(x = seq(xlimits[1], xlimits[2] - stratum_size_x, by=stratum_size_x),
                              y = seq(ylimits[1], ylimits[2] - stratum_size_y, by=stratum_size_y))

  generate_centers_within_stratum = function(min,max){
    runif(n = boxes_per_stratum, min, max)
  }

  # uniform random centers, with each stratum cell getting an equal number.
  center_x = purrr::map2(stratum_cells$x, stratum_cells$x+stratum_size_x, generate_centers_within_stratum)
  center_x = purrr::flatten_dbl(center_x)
  center_y = purrr::map2(stratum_cells$y, stratum_cells$y+stratum_size_y, generate_centers_within_stratum)
  center_y = purrr::flatten_dbl(center_y)

  boxes = dplyr::tibble(center_x = center_x,
                        center_y = center_y)

  total_boxes = nrow(boxes)

  boxes$box_id = 1:total_boxes
  boxes$size = box_size
  boxes$half_size = boxes$size/2

  return(boxes)
}


#' @description Fit the data to the described boxes using the estimator.
#'     This function does the spatial subsetting.
#'
#' @param data
#' @param boxes
#' @estimator function
#'
#' @return data.frame of length nrow(boxes), with estimates defined by
#'     estimator.
fit_estimators = function(data,
                          boxes,
                          estimator){

  run_estimator = function(box_id, size, half_size, center_x, center_y){
    data_subset = subset_points_to_box(data,
                                       box_center_x = center_x,
                                       box_center_y = center_y,
                                       box_half_size = half_size)

    subset_estimates = estimator(data_subset)
    subset_estimates$n_points = nrow(data_subset)
    subset_estimates$box_id = box_id
    return(subset_estimates)
  }

  # Get estimates of all metrics for each box
  box_estimates = purrr::pmap_df(boxes, run_estimator)

  fit_boxes = dplyr::left_join(boxes,box_estimates, by='box_id')

  return(fit_boxes)
}

#######################################
# Grid Helper Functions
#TODO: describe all these at once with @rdname
#######################################

#' @description \code{subset_points_to_box} subsets the points data.frame to those in the described box
#'
#' @param points data.frame
#'
#' @param box_center_x numeric
#'
#' @param box_center_y numeric
#'
#' @param box_half_size numeric
#'
#' @return data.frame subset of points
subset_points_to_box = function(points, box_center_x, box_center_y, box_half_size){
    dplyr::filter(points,
                   x > box_center_x - box_half_size,
                   x <= box_center_x + box_half_size,
                   y > box_center_y - box_half_size,
                   y <= box_center_y + box_half_size)
}

#' @description \code{subset_boxes_to_point} subsets boxes data.frame to those which contain the point x,y
#'
#' @param x numeric
#'
#' @param y numeric
#'
#' @param boxes data.frame
#'
#' @return data.frame subset of boxes
subset_boxes_to_point = function(x,y, boxes){
    dplyr::filter(boxes,
                 center_x - half_size < x,
                 center_x + half_size >= x,
                 center_y - half_size < y,
                 center_y + half_size >= y)
}

#' @description \code{within_bounds}: Is a number within some defined limits
#'
#' @param x numeric
#' @param low numeric
#' @param high numeric
#'
#' @return logical
within_bounds = function(x, low, high){
  x >= low & x <= high
}

#' @description \code{within_bounds2}: Is a number within some defined limits in 2 dimensions
#'
#' @param x numeric
#' @param x_low numeric
#' @param x_high numeric
#' @param y numeric
#' @param y_low numeric
#' @param y_high numeric
#'
#' @return logical
within_bounds2 = function(x, y,
                          x_low, x_high,
                          y_low, y_high){
  within_bounds(x, x_low, x_high) & within_bounds(y, y_low, y_high)
}
