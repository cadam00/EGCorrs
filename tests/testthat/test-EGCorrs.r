test_that("EGCorrs works", {

  library(EGCorrs)
  library(sf)
  set.seed(42, "Mersenne-Twister", sample.kind="Rejection")
  component_u       <- get_component_u()
  component_v       <- get_component_v()
  origin_areas      <- get_origin_areas()
  destination_areas <- get_destination_areas()
  out               <- suppressWarnings(
                           EGCorrs(component_u, component_v, origin_areas,
                           destination_areas, npoints = 1, lambda = 1,
                           niters = 10, k_neighbors = 7, nearest_grid_nodes = 4,
                           mask_shapefile = NULL, all_networks = FALSE)
  )

  expect_equal(length(out), 8)

  expect_equal(head(out$metrics_df$niters), c(1,2,3,4,5,6))

  expect_equal(out$list_of_networks, NULL)

  expect_equal(sf::st_bbox(out$solution_edges)[1:4],
               c(xmin = 22.87500397191,
                 ymin = 36.06250047306,
                 xmax = 24.45833743675,
                 ymax = 36.72916719341
               )
  )

  out               <- suppressWarnings(
    EGCorrs(component_u, component_v, origin_areas,
             destination_areas, npoints = 1, lambda = 1,
             niters = 10, k_neighbors = 7, nearest_grid_nodes = 4,
             mask_shapefile = NULL, all_networks = TRUE)
  )
  expect_equal(length(out$list_of_networks), 10)

  ## Examples of mask usage
  mask_shapefile <- st_as_sf(st_as_sfc(st_bbox(component_u),
                                       crs=terra::crs(component_u)))
  mask_shapefile <- st_crop(mask_shapefile, terra::ext(mask_shapefile) / 1.2)

  masked_result <- suppressWarnings(
                     EGCorrs(component_u, component_v, origin_areas,
                              destination_areas, npoints = 1, lambda = 1,
                              niters = 10, k_neighbors = 7,
                              nearest_grid_nodes = 4,
                              mask_shapefile = mask_shapefile,
                              all_networks = FALSE)
                   )

  expect_equal(sf::st_bbox(masked_result$solution_edges)[1:4],
    c(xmin = 23.00000398229,
      ymin = 36.10416714308,
      xmax = 24.33333742637,
      ymax = 36.68750052339
    )
  )

  # Check working under warnings
  vect_mask_shapefile <- terra::vect(mask_shapefile)
  masked_result <- suppressWarnings(
                      EGCorrs(component_u, component_v, origin_areas,
                               destination_areas, npoints = 1, lambda = 1,
                               niters = 10, k_neighbors = 7,
                               nearest_grid_nodes = 4,
                               mask_shapefile = vect_mask_shapefile,
                               all_networks = FALSE)
                    )

  expect_equal(sf::st_bbox(masked_result$solution_edges)[1:4],
               c(xmin = 23.00000398229,
                 ymin = 36.10416714308,
                 xmax = 24.33333742637,
                 ymax = 36.68750052339
               )
  )

  vect_mask_shapefile <- terra::project(vect_mask_shapefile, "+init=EPSG:4269")
  masked_result <- suppressWarnings(
    EGCorrs(component_u, component_v, origin_areas,
             destination_areas, npoints = 1, lambda = 1,
             niters = 10, k_neighbors = 7,
             nearest_grid_nodes = 4,
             mask_shapefile = vect_mask_shapefile,
             all_networks = FALSE)
  )
  expect_equal(sf::st_bbox(masked_result$solution_edges)[1:4],
               c(xmin = 23.00000398229,
                 ymin = 36.10416714308,
                 xmax = 24.33333742637,
                 ymax = 36.68750052339
               )
  )

  out               <- suppressWarnings(
    EGCorrs(component_u, component_v, origin_areas,
             destination_areas, npoints = 1, lambda = 1,
             niters = 10, k_neighbors = 7, nearest_grid_nodes = 4,
             mask_shapefile = NULL, all_networks = FALSE, progBar=FALSE)
  )


  ## Check warnings
  withWarnings <- function(expr) {
    myWarnings <- NULL
    wHandler <- function(w) {
      myWarnings <<- c(myWarnings, list(w))
      invokeRestart("muffleWarning")
    }
    val <- withCallingHandlers(expr, warning = wHandler)
    list(value = val, warnings = myWarnings)
  }

  masked_result <- suppressWarnings(
    withWarnings(EGCorrs(component_u, component_v,
                          st_transform(origin_areas,
                                       crs = "+init=EPSG:4269"),
                          st_transform(destination_areas,
                                       crs = "+init=EPSG:4269"),
                          npoints = 1, lambda = 1,
                          niters = 2, k_neighbors = 7,
                          nearest_grid_nodes = 4,
                          mask_shapefile =
                            st_transform(mask_shapefile,
                                         crs = "+init=EPSG:4269"),
                          all_networks = FALSE)))$warnings
  expect_equal(is(masked_result[[1]],"simpleWarning"), TRUE)
  expect_equal(is(masked_result[[2]],"simpleWarning"), TRUE)
  expect_equal(is(masked_result[[3]],"simpleWarning"), TRUE)
  expect_equal(is(masked_result[[4]],"simpleWarning"), TRUE)
  expect_equal(is(masked_result[[5]],"simpleWarning"), TRUE)
  expect_equal(masked_result[[1]]$message,
               paste0("Different resolution or coordinates among component_u",
              " and origin_areas:\nproject crs(component_u) on origin_areas."))
  expect_equal(masked_result[[2]]$message,
               paste0("Different resolution or coordinates among component_u",
    " and destination_areas:\nproject crs(component_u) on destination_areas."))
  expect_equal(masked_result[[3]]$message,
               paste0("Different resolution or coordinates among component_u",
          " and mask_shapefile:\nproject crs(component_u) on mask_shapefile"))
  expect_equal(masked_result[[4]]$message,
             paste0("attribute variables are assumed to be spatially constant",
                      " throughout all geometries"))
  expect_equal(masked_result[[5]]$message,
             paste0("attribute variables are assumed to be spatially constant",
                      " throughout all geometries"))

  out <- suppressWarnings(
    withWarnings(EGCorrs(component_u, component_v,
                          rbind(origin_areas, origin_areas),
                          rbind(destination_areas, destination_areas),
                          npoints = 1, lambda = 1,
                          niters = 2, k_neighbors = 7,
                          nearest_grid_nodes = 4,
                          mask_shapefile = NULL,
                          all_networks = FALSE)))$warnings
  expect_equal(is(out[[1]],"simpleWarning"), TRUE)
  expect_equal(out[[1]]$message,
               paste0("npoints < nrow(origin_areas).",
                      " Some origin areas will not have any point."))
  expect_equal(is(out[[2]],"simpleWarning"), TRUE)
  expect_equal(out[[2]]$message,
             paste0("npoints < nrow(destination_areas).",
                    " Some destination areas will not have any point."))

  ## Check errors
  masked_result <- try(suppressWarnings(
    EGCorrs(component_u, component_v, origin_areas,
                                destination_areas, npoints = 1, lambda = 1,
                                niters = 10, k_neighbors = 7,
                                nearest_grid_nodes = 4,
                                mask_shapefile = "Hello",
                                all_networks = FALSE)
                       ), silent = TRUE)

  expect_equal(class(masked_result) == "try-error", TRUE)

  out <- try(suppressWarnings(
     EGCorrs(component_u, component_v, origin_areas,
                                destination_areas, npoints = 1, lambda = 1,
                                niters = 10, k_neighbors = 7,
                                nearest_grid_nodes = "Hello",
                                mask_shapefile = NULL,
                                all_networks = FALSE)
                       ), silent = TRUE)

  expect_equal(class(out) == "try-error", TRUE)

  out <- try(suppressWarnings(
    EGCorrs(component_u, component_v, origin_areas,
                      destination_areas, npoints = 1, lambda = 1,
                      niters = "Hello", k_neighbors = 7,
                      nearest_grid_nodes = 4,
                      mask_shapefile = NULL,
                      all_networks = FALSE)
             ), silent = TRUE)

  expect_equal(class(out) == "try-error", TRUE)

  out <- try(suppressWarnings(
    EGCorrs(component_u, component_v, origin_areas,
                      destination_areas, npoints = 1, lambda = 1,
                      niters = 10, k_neighbors = "Hello",
                      nearest_grid_nodes = 4,
                      mask_shapefile = NULL,
                      all_networks = FALSE)
             ), silent = TRUE)

  expect_equal(class(out) == "try-error", TRUE)

  out <- try(suppressWarnings(
    EGCorrs(component_u="Hello", component_v, origin_areas,
             destination_areas, npoints = 1, lambda = 1,
             niters = 10, k_neighbors = 7,
             nearest_grid_nodes = 4,
             mask_shapefile = NULL,
             all_networks = FALSE)
  ), silent = TRUE)

  expect_equal(class(out) == "try-error", TRUE)

  out <- try(suppressWarnings(
    EGCorrs(component_u, component_v="Hello", origin_areas,
             destination_areas, npoints = 1, lambda = 1,
             niters = 10, k_neighbors = 7,
             nearest_grid_nodes = 4,
             mask_shapefile = NULL,
             all_networks = FALSE)
  ), silent = TRUE)

  expect_equal(class(out) == "try-error", TRUE)

  out <- try(suppressWarnings(
    EGCorrs(component_u, component_v,
             origin_areas=matrix(c("Hello","Hello"), nrow=2),
             destination_areas, npoints = 1, lambda = 1,
             niters = 10, k_neighbors = 7,
             nearest_grid_nodes = 4,
             mask_shapefile = NULL,
             all_networks = FALSE)
  ), silent = TRUE)

  expect_equal(class(out) == "try-error", TRUE)

  out <- try(suppressWarnings(
    EGCorrs(component_u, component_v, origin_areas,
             destination_areas=matrix(c("Hello","Hello"), nrow=2),
             npoints = 1, lambda = 1,
             niters = 10, k_neighbors = 7,
             nearest_grid_nodes = 4,
             mask_shapefile = NULL,
             all_networks = FALSE)
  ), silent = TRUE)

  expect_equal(class(out) == "try-error", TRUE)

  out <- try(suppressWarnings(
    EGCorrs(component_u, component_v, origin_areas,
                      destination_areas, npoints = 1.5, lambda = 1,
                      niters = 10, k_neighbors = 7,
                      nearest_grid_nodes = 4,
                      mask_shapefile = NULL,
                      all_networks = FALSE)
             ), silent = TRUE)

  expect_equal(class(out) == "try-error", TRUE)

  out <- try(suppressWarnings(
    EGCorrs(component_u, component_v, origin_areas,
                      destination_areas, npoints = 1, lambda = 1,
                      niters = 10, k_neighbors = 7.5,
                      nearest_grid_nodes = 4,
                      mask_shapefile = NULL,
                      all_networks = FALSE)
             ), silent = TRUE)

  expect_equal(class(out) == "try-error", TRUE)

  out <- try(suppressWarnings(
    EGCorrs(component_u, component_v, origin_areas,
                      destination_areas, npoints = 1, lambda = 1,
                      niters = 10.5, k_neighbors = 7,
                      nearest_grid_nodes = 4,
                      mask_shapefile = NULL,
                      all_networks = FALSE)
             ), silent = TRUE)

  expect_equal(class(out) == "try-error", TRUE)

  out <- try(suppressWarnings(
    EGCorrs(component_u, component_v, origin_areas,
                      destination_areas, npoints = 1, lambda = 1,
                      niters = 10, k_neighbors = 7,
                      nearest_grid_nodes = 4.5,
                      mask_shapefile = NULL,
                      all_networks = FALSE)
             ), silent = TRUE)

  expect_equal(class(out) == "try-error", TRUE)

  out <- try(suppressWarnings(
    EGCorrs(component_u, component_v, origin_areas,
                      destination_areas, npoints = -1, lambda = 1,
                      niters = 10, k_neighbors = 7,
                      nearest_grid_nodes = 4,
                      mask_shapefile = NULL,
                      all_networks = FALSE)
             ), silent = TRUE)

  expect_equal(class(out) == "try-error", TRUE)

  out <- try(suppressWarnings(
    EGCorrs(component_u, component_v, origin_areas,
                      destination_areas, npoints = 1, lambda = 1,
                      niters = 10, k_neighbors = -7,
                      nearest_grid_nodes = 4,
                      mask_shapefile = NULL,
                      all_networks = FALSE)
             ), silent = TRUE)

  expect_equal(class(out) == "try-error", TRUE)

  out <- try(suppressWarnings(
    EGCorrs(component_u, component_v, origin_areas,
                      destination_areas, npoints = 1, lambda = 1,
                      niters = -10, k_neighbors = 7,
                      nearest_grid_nodes = 4,
                      mask_shapefile = NULL,
                      all_networks = FALSE)
             ), silent = TRUE)

  expect_equal(class(out) == "try-error", TRUE)

  out <- try(suppressWarnings(
    EGCorrs(component_u, component_v, origin_areas,
                      destination_areas, npoints = 1, lambda = 1,
                      niters = 10, k_neighbors = 7,
                      nearest_grid_nodes = -4,
                      mask_shapefile = NULL,
                      all_networks = FALSE)
             ), silent = TRUE)

  expect_equal(class(out) == "try-error", TRUE)

  out <- try(suppressWarnings(
    EGCorrs(component_u, component_v, st_cast(origin_areas,"LINESTRING"),
             destination_areas, npoints = 1, lambda = 1,
             niters = 10, k_neighbors = 7,
             nearest_grid_nodes = 4,
             mask_shapefile = NULL,
             all_networks = FALSE)
  ), silent = TRUE)

  expect_equal(class(out) == "try-error", TRUE)

  out <- try(suppressWarnings(
    EGCorrs(component_u, component_v, origin_areas,
             st_cast(destination_areas,"LINESTRING"), npoints = 1, lambda = 1,
             niters = 10, k_neighbors = 7,
             nearest_grid_nodes = 4,
             mask_shapefile = NULL,
             all_networks = FALSE)
  ), silent = TRUE)

  expect_equal(class(out) == "try-error", TRUE)

  terra::crs(component_v) <- ""
  out <- try(suppressWarnings(
    EGCorrs(component_u, component_v, origin_areas,
                      destination_areas, npoints = 1, lambda = 1,
                      niters = 10, k_neighbors = 7,
                      nearest_grid_nodes = 4,
                      mask_shapefile = NULL,
                      all_networks = FALSE)
             ), silent = TRUE)

  expect_equal(class(out) == "try-error", TRUE)

  component_v <- get_component_v()
  terra::res(component_v) <- c(1, 1)
  terra::ext(component_v) <- terra::ext(component_u)
  out <- try(suppressWarnings(
    EGCorrs(component_u, component_v, origin_areas,
                      destination_areas, npoints = 1, lambda = 1,
                      niters = 10, k_neighbors = 7,
                      nearest_grid_nodes = 4,
                      mask_shapefile = NULL,
                      all_networks = FALSE)
             ), silent = TRUE)

  expect_equal(class(out) == "try-error", TRUE)

  component_v <- get_component_v()
  component_v <- terra::crop(component_v, terra::ext(component_v) / 4)
  out <- try(suppressWarnings(
    EGCorrs(component_u, component_v, origin_areas,
                      destination_areas, npoints = 1, lambda = 1,
                      niters = 10, k_neighbors = 7,
                      nearest_grid_nodes = 4,
                      mask_shapefile = NULL,
                      all_networks = FALSE)
             ), silent = TRUE)

  expect_equal(class(out) == "try-error", TRUE)

})
