context('Testing the flowering gradient simulator spatial_flowering_sampler()')

x = seq(0,1,0.05)
y = seq(0,1,0.05)

test_that('random gradient generators are reproducible given a seed', {
  run_1 = generate_random_spatial_gradient(x,y,seed=5)
  expect_equal(generate_random_spatial_gradient(x,y,seed=5), run_1)
  
  run_2 = generate_random_spatial_gradient(x,y,seed=6)
  expect_equal(generate_random_spatial_gradient(x,y,seed=6), run_2)
  
  run_3 = generate_random_spatial_gradient(x,y,seed=6000)
  expect_equal(generate_random_spatial_gradient(x,y,seed=6000), run_3)
  
  run_4 = spatial_flowering_grid(xlimits=c(0,1),
                                 ylimits = c(0,1),
                                 cell_size=0.4,
                                 start_doy=90,
                                 flowering_length = 30,
                                 flowering_gradient = 10/0.1,
                                 spatial_gradient_type = 'non-linear',
                                 seed=15)
  run_5 = spatial_flowering_grid(xlimits=c(0,1),
                                 ylimits = c(0,1),
                                 cell_size=0.4,
                                 start_doy=90,
                                 flowering_length = 30,
                                 flowering_gradient = 10/0.1,
                                 spatial_gradient_type = 'non-linear',
                                 seed=15)
  expect_equal(run_4, run_5)
})

test_that('primary simulation function has correct output',{
  sample_size = 100
  simulated_data = spatial_flowering_sampler(n=sample_size)
  
  # sample size specified is correct
  expect_equal(nrow(simulated_data), sample_size)
  # no na's
  expect_true(all(!is.na(simulated_data)))
  
  # should be flowering presence only
  simulated_data_po = spatial_flowering_sampler(n=sample_size, sample_type = 'po')
  expect_true(all(simulated_data_po$flower_present==1))
  
  # should have flowering presence and absence, and in the proportion specified
  fraction_present = 0.5
  simulated_data_pa = spatial_flowering_sampler(n=sample_size, 
                                                sample_type = 'pa',
                                                fraction_present = fraction_present)
  expect_equal(sort(unique(simulated_data_pa$flower_present)), c(0,1))
  expect_equal(mean(simulated_data_pa$flower_present), fraction_present)
  
  # DOY is always between 1-365
  expect_true(all(simulated_data$doy >= 1))
  expect_true(all(simulated_data$doy <= 365))
  expect_true(all(simulated_data_po$doy >= 1))
  expect_true(all(simulated_data_po$doy <= 365))
  # Need to fix this in pa generation, see https://github.com/sdtaylor/phenology_gradients/issues/3
  #expect_true(all(simulated_data_pa$doy >= 1))
  #expect_true(all(simulated_data_pa$doy <= 365))
})
