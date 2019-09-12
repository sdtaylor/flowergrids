context('Testing the spatial weibull model weibull_grid()')

test_that('model function runs', {
  sample_points = spatial_flowering_sampler(n=200,
                                          xlimits = c(0,1),
                                          ylimits = c(0,1),
                                          start_doy = 90,
                                          flowering_length = 30,
                                          flowering_gradient = 10/0.1,
                                          spatial_gradient_type = 'non-linear',
                                          clustering = FALSE,
                                          sample_type = 'po',
                                          seed=100)
  
  fitted_model = weibull_grid(doy_points = sample_points,
                                      stratum_size_x = 0.1,
                                      stratum_size_y = 0.1,
                                      boxes_per_stratum = 5,
                                      box_size = 0.2,
                                      xlimits = c(0,1),
                                      ylimits = c(0,1),
                                      edge_buffer = 0)  
  
  new_points = expand.grid(x=seq(0,1,0.2), y=seq(0,1,0.2))
  
  predictions = predict.weibull_grid(fitted_model, new_points, type='onset')
  
  expect_equal(length(predictions), nrow(new_points))
})

test_that('model helper function subset_points_to_boxes', {
  test_points = round(expand.grid(x=seq(0,1,0.2),y=seq(0,1,0.2)),1)
  
  subset_points_1 = subset_points_to_box(test_points, box_center_x = 0.5, box_center_y = 0.25, box_half_size = 0.2)
  expect_true(all(subset_points_1 == dplyr::tibble(x=c(0.4,0.6,0.4,0.6), y=c(0.2,0.2,0.4,0.4))))
  
  subset_points_2 = subset_points_to_box(test_points, box_center_x = 0.3, box_center_y = 0.7, box_half_size = 0.12)
  expect_true(all(subset_points_2 == dplyr::tibble(x=c(0.2,0.4,0.2,0.4), y=c(0.6,0.6,0.8,0.8))))
  
  subset_points_3 = subset_points_to_box(test_points, box_center_x = 0.8, box_center_y = 0.8, box_half_size = 0.05)
  expect_true(all(subset_points_3 == dplyr::tibble(x=c(0.8), y=c(0.8))))
})

test_that('model helper function subset_boxes_to_points', {
  test_boxes = dplyr::tibble(center_x=c(0.1, 0.3, 0.5, 0.5), 
                             center_y = c(0.7, 0.5, 0.3, 0.8), 
                             half_size=c(0.3, 0.2, 0.1, 0.2),
                             box_id = c(1,2,3,4))
  
  subset_boxes_1 = subset_boxes_to_point(x = 0.2, y=0.6, test_boxes)
  expect_equal(subset_boxes_1$box_id, c(1,2))
  
  subset_boxes_2 = subset_boxes_to_point(x = 0.55, y=0.3, test_boxes)
  expect_equal(subset_boxes_2$box_id, c(3))
  
  subset_boxes_3 = subset_boxes_to_point(x = 0.35, y=0.65, test_boxes)
  expect_equal(subset_boxes_3$box_id, c(1,2,4))
  
  # Should be no boxes matching this point
  subset_boxes_4 = subset_boxes_to_point(x = 0.6, y=0.5, test_boxes)
  expect_equal(nrow(subset_boxes_4), 0)
})

test_that('bounding functions work', {
  expect_true(within_bounds(10,5,15))
  expect_true(within_bounds(0.1,0,1))
  expect_true(within_bounds(0.5,0.5,0.51))
  expect_false(within_bounds(10,15,20))
  expect_false(within_bounds(1,1.1,20))
  expect_false(within_bounds(0.5,0.51,1))
  
  expect_true(within_bounds2(x=10, y=5,
                             x_low=5,x_high = 15,
                             y_low=1,y_high = 10))
  expect_true(within_bounds2(x=0.1, y=0.8,
                             x_low=0,x_high = 1,
                             y_low=0.5,y_high = 1))
  expect_true(within_bounds2(x=0.1, y=0.1,
                             x_low=0,x_high = 0.1,
                             y_low=0,y_high = 0.11))
  expect_false(within_bounds2(x=1, y=1,
                             x_low=0,x_high = 0.9,
                             y_low=0,y_high = 0.9))
  expect_false(within_bounds2(x=0.5, y=0.5,
                               x_low=0,x_high = 1,
                               y_low=0,y_high = 0.49))
})  

