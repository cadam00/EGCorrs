test_that("get_origin_areas works", {
  origin_areas <- EGCorrs::get_origin_areas()
  expect_equal(sf::st_bbox(origin_areas)[1:4],
               c(xmin = 22.86968610966647,
                 ymin = 36.04912309180337,
                 xmax = 23.13855293332797,
                 ymax = 36.39421267146248
               )
  )
  expect_equal(dim(origin_areas), c(1, 2))
})


