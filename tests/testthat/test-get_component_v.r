test_that("get_component_v works", {
  component_v <- EGCorrs::get_component_u()
  expect_equal(terra::ext(component_v)[1:4],
               c(xmin = 22.8541706368489,
                 xmax = 24.4791707718168,
                 ymin = 36.0416671380519,
                 ymax = 36.7500005284211
               )
  )
  expect_equal(dim(component_v), c(17, 39, 1))
})
