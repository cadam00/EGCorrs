test_that("get_destination_areas works", {
  destination_areas <- EGCorrs::get_destination_areas()
  expect_equal(sf::st_bbox(destination_areas)[1:4],
               c(xmin = 24.30529576997133,
                 ymin = 36.62172860198531,
                 xmax = 24.49323320551187,
                 ymax = 36.75966068439976
               )
  )
  expect_equal(dim(destination_areas), c(1, 2))
})


