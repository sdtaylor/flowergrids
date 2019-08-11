#' @description The Weibull Grid method as described in Taylor et al.
#'       Combines the methodology of Fink et al. 2010 with the Weibull estimator
#'       as described in Pearse et al. 2017
#'

weibull_grid = function(doy_points,
                        stratum_size_x=0.1,
                        stratum_size_y=0.1,
                        boxes_per_stratum=5,
                        box_size=0.3,
                        xlimits=c(0,1),
                        ylimits=c(0,1),
                        edge_buffer=0.1,
                        not_enough_data_fallback='use_na',
                        max_n_per_box=50){

  model_details = list()
  model_details$doy_points = doy_points
  model_details$stratum_size_x = stratum_size_x
  model_details$stratum_size_y = stratum_size_y
  model_details$boxes_per_stratum = boxes_per_stratum
  model_details$box_size = box_size
  model_details$xlimits = xlimits
  model_details$ylimits = ylimits
  model_details$edge_buffer = edge_buffer
  model_details$max_n_per_box = max_n_per_box
  model_details$not_enough_data_fallback = not_enough_data_fallback

  # The uniformly random boxes.
  boxes = create_grid_boxes(stratum_size_x = stratum_size_x,
                            stratum_size_y = stratum_size_y,
                            boxes_per_stratum = boxes_per_stratum,
                            box_size = box_size,
                            xlimits = xlimits,
                            ylimits = ylimits)

  # This take a subset of doy_points and return the estimated
  # values
  weibull_estimator_for_grid = function(doy_points_subset){
    estimates = list()
    # estimates$onset_estimate = tryCatch(as.numeric(phest::weib.limit(doy_points_subset$doy, k=max_n_per_box)[1]),
    #                                     error = function(cond){fallback(doy_points_subset$doy)})
    # estimates$end_estimate   = tryCatch(as.numeric(phest::weib.limit(doy_points_subset$doy*-1, k=max_n_per_box)[1]) * -1,
    #                                     error = function(cond){fallback(doy_points_subset$doy * -1) * -1})
    estimates$onset_estimate = phest::weib.limit(doy_points_subset$doy, k=50)[1]
    estimates$end_estimate   = as.numeric(phest::weib.limit(doy_points_subset$doy*-1, k=50)[1]) * -1
    estimates$peak_estimate  = mean(doy_points_subset$doy)
    return(estimates)
  }

  # fit_estimators does the spatial subsetting and fitting
  model_details$fitted_boxes = fit_estimators(boxes = boxes,
                                              data  = doy_points,
                                              estimator = weibull_estimator_for_grid)

  return(structure(model_details, class = 'weibull_grid'))
}

#' @description Make predictions using a Weibull Grid model
#'
#' @param model
#' @param doy_points data.frame
#' @param type
#' @param se
#' @param level
#'
#' @return vector if se is FALSE, data.frame if TRUE
predict.weibull_grid = function(model,
                                doy_points,
                                type = 'onset',
                                se = F,
                                level = 0.95){

  outside_buffer = function(x,y){
    !within_bounds2(x,y,
                    x_low  = model$xlimits[1] + model$edge_buffer,
                    x_high = model$xlimits[2] - model$edge_buffer,
                    y_low  = model$ylimits[1] + model$edge_buffer,
                    y_high = model$ylimits[2] - model$edge_buffer)
  }

  lower_quantile = (1 - level)/2
  upper_quantile = 1 - lower_quantile

  estimate_metrics_from_model = function(x, y){
    box_subset = subset_boxes_to_point(x = x,
                                       y = y,
                                       boxes = model$fitted_boxes)

    estimates = list()
    estimates$x = x
    estimates$y = y
    if(outside_buffer(x,y)){
      estimates$onset_estimate       = NA
      estimates$onset_estimate_upper = NA
      estimates$onset_estimate_lower = NA

      estimates$end_estimate         = NA
      estimates$end_estimate_upper   = NA
      estimates$end_estimate_lower   = NA

      estimates$peak_estimate        = NA
      estimates$peak_estimate_upper  = NA
      estimates$peak_estimate_lower  = NA

      estimates$outside_buffer       = TRUE
    } else {
      estimates$onset_estimate       = median(box_subset$onset_estimate, na.rm=T)
      estimates$onset_estimate_upper = quantile(box_subset$onset_estimate, upper_quantile, na.rm=T)
      estimates$onset_estimate_lower = quantile(box_subset$onset_estimate, lower_quantile, na.rm=T)

      estimates$end_estimate         = median(box_subset$end_estimate, na.rm=T)
      estimates$end_estimate_upper   = quantile(box_subset$end_estimate, upper_quantile, na.rm=T)
      estimates$end_estimate_lower   = quantile(box_subset$end_estimate, lower_quantile, na.rm=T)

      estimates$peak_estimate        = median(box_subset$peak_estimate, na.rm=T)
      estimates$peak_estimate_upper  = quantile(box_subset$peak_estimate, upper_quantile, na.rm=T)
      estimates$peak_estimate_lower  = quantile(box_subset$peak_estimate, lower_quantile, na.rm=T)

      estimates$outside_buffer       = FALSE
    }
    return(estimates)
  }

  # Get estimates for each prediction point
  point_estimates = purrr::pmap_df(doy_points[c('x','y')], estimate_metrics_from_model)

  outside_buffer_count = sum(point_estimates$outside_buffer)
  if(outside_buffer_count>0){
    warning(paste(outside_buffer_count,'points were outside the buffer and could not be estimated.'))
  }

  if(type == 'onset'){
    point_estimates$estimate       = point_estimates$onset_estimate
    point_estimates$estimate_lower = point_estimates$onset_estimate_lower
    point_estimates$estimate_upper = point_estimates$onset_estimate_upper
  } else if(type == 'end'){
    point_estimates$estimate       = point_estimates$end_estimate
    point_estimates$estimate_lower = point_estimates$end_estimate_lower
    point_estimates$estimate_upper = point_estimates$end_estimate_upper
  } else if(type == 'peak'){
    point_estimates$estimate       = point_estimates$peak_estimate
    point_estimates$estimate_lower = point_estimates$peak_estimate_lower
    point_estimates$estimate_upper = point_estimates$peak_estimate_upper
  } else {
    stop(paste('unknown prediction type: ',type))
  }

  if(se){
    return(dplyr::select(point_estimates, estimate, estimate_lower, estimate_upper))
  } else{
    return(point_estimates$estimate)
  }
}
