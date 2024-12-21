get_component_u <- function()
  rast(system.file("external/component_u.tiff", package="EGCorrs"))
get_component_v <- function()
  rast(system.file("external/component_v.tiff", package="EGCorrs"))
get_origin_areas <- function(){
  origin_areas <-
    st_read(system.file("external/origin_areas/origin_areas.shp",
                      package="EGCorrs"))
  st_transform(origin_areas, crs(get_component_u()))
}

get_destination_areas <- function(){
  destination_areas <-
    st_read(system.file("external/destination_areas/destination_areas.shp",
                        package="EGCorrs"))
  st_transform(destination_areas, crs(get_component_u()))
}

points_function<-function(npoints, areas){
  # Create random points inside these multipolygons
  n_multipolygons <- min(nrow(areas), npoints)
  outpoints <- vector("list", n_multipolygons)

  points_vec <- rep( npoints / n_multipolygons, n_multipolygons)
  points_new <- points_vec
  points_new <- floor(points_vec)
  indices    <- tail(order(points_vec-points_new),
                     round(sum(points_vec)) - sum(points_new))
  points_new[indices] <- points_new[indices] + 1

  for (i in seq_len(n_multipolygons)){
    outpoints[[i]] <-
      st_sample(st_make_valid(areas[i,]),
                size = points_new[i],
                type = "random")
  }
  return(do.call(c, outpoints))
}
melt_neighbours <- function(input_to_melt){
  dim_input_to_melt <- dim(input_to_melt)
  nn_results_edges <- matrix(c(rep(seq_len(dim_input_to_melt[1]),
                                   times = dim_input_to_melt[2]),
                               input_to_melt),
                             nrow = dim_input_to_melt[1] * dim_input_to_melt[2],
                             ncol = 2
  )
  colnames(nn_results_edges) <- c("from", "to")
  nn_results_edges <-
    nn_results_edges[nn_results_edges[,1] != nn_results_edges[,2],]
  return(cbind(nn_results_edges[!duplicated(
    t(vapply(seq_len(nrow(nn_results_edges)),
             function(i)
               sort(nn_results_edges[i,]),
             numeric(2)))),],
    "weight" = 1)
  )
}
# Function to normalize a column to the range [0.0..., 1]
normalize_dividing <- function(x) {
  #(x - min(x)) / (max(x) - min(x))
  #(x - 0) / (max(x) - 0) #we consider as the minimum value being 0
  # and all of them being positive integers
  x/max(x)
}
greedy_match<- function(dist_matrix) {

  dim_dist        <- dim(dist_matrix)
  rownames_dist   <- rownames(dist_matrix)
  colnames_dist   <- colnames(dist_matrix)
  matched_pairs   <- matrix(nrow = dim_dist[1], ncol = 2)
  col_unavailable <- logical(dim_dist[2])

  for (cc in seq_len(dim_dist[1])) {

    # From cc point to all the ending points
    current_row                  <- dist_matrix[cc, ]

    # Set unavailable columns to Inf in the current row
    current_row[col_unavailable] <- Inf

    # Find the column with the minimum distance
    min_dist_col <- which.min(current_row)

    # Store the matched pair
    matched_pairs[cc, ] <- c(rownames_dist[cc], colnames_dist[min_dist_col])

    # Mark the column as unavailable
    col_unavailable[min_dist_col] <- TRUE

  }

  return(matched_pairs)
}

EGCorrs <- function(component_u, component_v, origin_areas, destination_areas,
                     npoints = 40, lambda = 1, niters = 100, k_neighbors = 7,
                     nearest_grid_nodes = 4, mask_shapefile = NULL,
                     all_networks = FALSE, progBar = TRUE){

  # Several parts of the code for Sea Currents to Connectivity Transformation
  # are taken from the SeaGraphs R package
  # https://github.com/cadam00/SeaGraphs

  if (!is(component_u, "SpatRaster"))
    stop("component_u must be a SpatRaster object.")
  if (!is(component_v, "SpatRaster"))
    stop("component_v must be a SpatRaster object.")

  if ( (npoints <= 0) )
    stop("npoints must must be a positive number.")
  if (npoints != as.integer(npoints))
    stop("npoints must be an integer number.")
  if ( npoints <  nrow(origin_areas))
    warning(
      paste0("npoints < nrow(origin_areas). Some origin areas will",
             " not have any point.")
    )
  if ( npoints <  nrow(destination_areas))
    warning(
      paste0("npoints < nrow(destination_areas). Some destination areas will",
             " not have any point.")
    )

  if ( (k_neighbors <= 0) )
    stop("k_neighbors must must be a positive number.")
  if (k_neighbors != as.integer(k_neighbors))
    stop("k_neighbors must be an integer number.")

  if ( (nearest_grid_nodes <= 0) )
    stop("nearest_grid_nodes must must be a positive number.")
  if (nearest_grid_nodes != as.integer(nearest_grid_nodes))
    stop("nearest_grid_nodes must be an integer number.")

  if ( (niters < 0) )
    stop("niters must must be a non-negative number.")
  if (niters != as.integer(niters))
    stop("niters must be an integer number.")

  crs_component_u <- crs(component_u)
  crs_component_v <- crs(component_v)

  if (crs_component_u != crs_component_v)
    stop(
      paste0("Different coordinate system between component_u and component_v.",
             " Check with crs(component_u) and crs(component_v).")

    )

  if (ext(component_u) != ext(component_v))
    stop(
      paste0("Different extent between component_u and component_v.",
             " Check with ext(component_u) and ext(component_v).")

    )

  if (any(res(component_u) != res(component_v)))
    stop(
      paste0("Different resolution between component_u and component_v.",
             " Check with res(component_u) and res(component_v).")

    )

  if (!((is(origin_areas, "sf") && is(origin_areas, "data.frame")) |
        is(origin_areas, "sfc")))
    stop(
      paste0("origin_areas arguement must be both sf and data.frame or sfc.")

    )
  if (!((is(destination_areas, "sf") && is(destination_areas, "data.frame")) |
        is(destination_areas, "sfc")))
    stop(
    paste0("destination_areas arguement must be both sf and data.frame or sfc.")

    )

  if (!all(st_is(origin_areas, 'MULTIPOLYGON') |
           st_is(origin_areas, 'sfc_MULTIPOLYGON') |
           st_is(origin_areas, 'POLYGON') |
           st_is(origin_areas, 'sfc_POLYGON')
      ))
    stop(
      "origin_areas must be rows of 'POLYGON' or'MULTIPOLYGON' geometry type.")
  if (!all(st_is(destination_areas, 'MULTIPOLYGON') |
           st_is(destination_areas, 'sfc_MULTIPOLYGON') |
           st_is(destination_areas, 'POLYGON') |
           st_is(destination_areas, 'sfc_POLYGON')
      ))
    stop(
  "destination_areas must be rows of 'POLYGON' or'MULTIPOLYGON' geometry type.")

  if (st_crs(origin_areas)$wkt != crs_component_u){
    warning(
      paste0(
        "Different resolution or coordinates among component_u and",
        " origin_areas:\nproject crs(component_u) on origin_areas."
      )
    )
    origin_areas <- st_transform(origin_areas, crs = crs_component_u)
  }
  if (st_crs(destination_areas)$wkt != crs_component_u){
    warning(
      paste0(
        "Different resolution or coordinates among component_u and",
        " destination_areas:\nproject crs(component_u) on destination_areas."
      )
    )
    destination_areas <- st_transform(destination_areas, crs = crs_component_u)
  }


  if (!is.null(mask_shapefile)){

    is_mask_shapefile_sf         <- is(mask_shapefile, "sf")
    is_mask_shapefile_SpatVector <- is(mask_shapefile, "SpatVector")

    if (is_mask_shapefile_sf) {
      if (st_crs(mask_shapefile)$wkt != crs_component_u){
        warning(
          paste0(
            "Different resolution or coordinates among component_u and",
            " mask_shapefile:\nproject crs(component_u) on mask_shapefile"
          )
        )
        mask_shapefile <- st_transform(mask_shapefile, crs = crs_component_u)
      }
    } else if (is_mask_shapefile_SpatVector) {
      if (crs(mask_shapefile) != crs_component_u){
        warning(
          paste0(
            "Different resolution or coordinates among component_u and",
            " mask_shapefile:\nproject crs(component_v) on mask_shapefile")
        )
        mask_shapefile <- project(mask_shapefile, crs_component_u)
      }

    } else {
      stop("mask_shapefile must be either 'sf' or 'SpatVector' class.")
    }
    component_u <- crop(component_u, mask_shapefile, mask=TRUE)
    component_v <- crop(component_v, mask_shapefile, mask=TRUE)
  }

  sf_u              <- st_as_sf(as.polygons(component_u))
  origin_areas      <- st_intersection(origin_areas,      sf_u)
  destination_areas <- st_intersection(destination_areas, sf_u)

  random_points_1 <- points_function(npoints, origin_areas)

  random_points_2 <- points_function(npoints, destination_areas)

  centroids_within_polygon <- st_as_sf(as.data.frame(crds(component_u)),
                                       coords = c("x","y"),
                                       crs=crs_component_u)[[1]]

  vect_centroids_within_polygon <- vect(c(centroids_within_polygon,
                                          random_points_1,
                                          random_points_2))
  points_component_u <- extract(component_u, vect_centroids_within_polygon)
  points_component_v <- extract(component_v, vect_centroids_within_polygon)


  two_ngnp1           <- seq.int(from=2, to=(nearest_grid_nodes+1), by=1)

  distances_1          <- st_distance(random_points_1, centroids_within_polygon)
  closest_metrics_dframe_1 <- matrix(ncol=nearest_grid_nodes,
                                     nrow=length(random_points_1))
  for (i in seq_along(random_points_1)){
    closest_metrics_dframe_1[i,] <- order(distances_1[i,])[two_ngnp1]
  }

  length_centroids_within_polygon <- length(centroids_within_polygon)

  edge_list_1     <- melt_neighbours(closest_metrics_dframe_1)
  edge_list_1[,1] <- edge_list_1[,1] + length_centroids_within_polygon

  distances_2          <- st_distance(random_points_2, centroids_within_polygon)
  closest_metrics_dframe_2 <- matrix(ncol=nearest_grid_nodes,
                                     nrow=length(random_points_2))
  for (i in seq_along(random_points_2)){
    closest_metrics_dframe_2[i,] <- order(distances_2[i,])[two_ngnp1]
  }

  edge_list_2     <- melt_neighbours(closest_metrics_dframe_2)
  edge_list_2[,1] <- edge_list_2[,1] + max(edge_list_1[,1])

  distances_nn      <- st_distance(centroids_within_polygon)
  nrow_distances_nn <- nrow(distances_nn)
  input_grid_graph  <- matrix(ncol = k_neighbors, nrow = nrow_distances_nn)
  two_kp1           <- seq.int(from=2, to=(k_neighbors+1), by=1)
  for (i in seq_len(nrow_distances_nn)){
    input_grid_graph[i,] <- order(distances_nn[i,])[two_kp1]
  }

  edge_list_0 <- melt_neighbours(input_grid_graph)

  #unifies edge list
  edge_list_i <- rbind(edge_list_0, edge_list_1, edge_list_2)

  points_object <- c(centroids_within_polygon, random_points_1, random_points_2)

  nrow_edge_list_i <- nrow(edge_list_i)
  edges_list       <- vector(mode = "list", length  = nrow_edge_list_i)
  values           <- numeric(nrow_edge_list_i)
  for (i in seq_len(nrow_edge_list_i)){
    first_point     <- edge_list_i[i,1]
    second_point    <- edge_list_i[i,2]
    edges_list[[i]] <- st_linestring(
      c(points_object[[first_point]],
        points_object[[second_point]])
    )

    xy_current_vector_first  <- c(points_component_u[first_point,2],
                                  points_component_v[first_point,2])
    xy_current_vector_second <- c(points_component_u[second_point,2],
                                  points_component_v[second_point,2])
    xy_current_vector        <- (xy_current_vector_first +
                                   xy_current_vector_second) / 2

    coords_basic_linestring <- st_coordinates(edges_list[[i]])

    # Calculate the direction vector
    xy_direction_vector <- coords_basic_linestring[2, c(1,2)] -
      coords_basic_linestring[1, c(1,2)]

    ## https://en.wikipedia.org/wiki/Vector_projection
    ##                                        #Definitions_in_terms_of_a_and_b
    ## Scalar projection of a_ onto b_ = ||a_|| * cos(theta) =
    ##                                                    dot(a_, b_) / ||b_||
    values[i] <- sum(xy_current_vector * xy_direction_vector) /
      sqrt(sum(xy_direction_vector ^ 2))
  }

  sfc_lines <- do.call(st_sfc, edges_list)

  values[is.na(values)]                  <-  0
  values[values >= 0  & values <  0.001] <-  0.001
  values[values <  0  & values > -0.001] <- -0.001

  sfc_lines        <- st_as_sf(sfc_lines)
  sfc_lines$weight <- values

  negative             <- sfc_lines$weight < 0
  sfc_lines[negative,] <- st_reverse(sfc_lines[negative,])

  sfc_lines$weight <- normalize_dividing(st_length(sfc_lines) /
                                           abs(sfc_lines$weight))
  sfc_lines$weight <- sfc_lines$weight * lambda

  net_result         <- as_sfnetwork(sfc_lines, directed=TRUE)
  st_crs(net_result) <- crs_component_u

  result_edges <- st_as_sf(net_result, "edges")
  result_nodes <- st_as_sf(activate(net_result, "nodes"))

  starting_points_ID <- seq.int(
    from = length_centroids_within_polygon + 1,
    to   = length_centroids_within_polygon + npoints
  )
  ending_points_ID   <- seq.int(
    from = length_centroids_within_polygon + npoints + 1,
    to   = length_centroids_within_polygon + (2 * npoints)
  )

  # starting_points <- result_nodes[starting_points_ID,]
  # ending_points   <- result_nodes[ending_points_ID,]

  initial_w <- edge_attr(net_result, "weight")
  initial_w_matrix <- data.frame(
    "link"   = seq_along(initial_w),
    "weight" = initial_w
  )

  population <- niters * npoints

  net_result_loop  <- net_result

  ###### LOOP!!!!
  list_of_edges_freq <- vector(mode = "list", length = niters)
  list_of_congestion <- vector(mode = "list", length = niters)

  if (all_networks){
    list_of_networks <- vector(mode = "list", length = niters)
  } else {
    list_of_networks <- NULL
  }

  if (progBar){
    pb <- txtProgressBar(0, length(niters), style = 3)
    on.exit(close(pb))
  }
  seq_int_niters  <- seq_len(niters)
  seq_len_npoints <- seq_len(npoints)
  for(i in seq_int_niters){
    if (progBar) setTxtProgressBar(pb, i)
    # if (verbose) message(paste0("Iteration: ", i))

    # Get the shortest distances from all the starting points to all the
    # ending points
    find_path_distances <- distances(graph = net_result_loop,
                                     v     = starting_points_ID,
                                     to    = ending_points_ID,
                                     mode  = "out"
    )

    colnames(find_path_distances) <- ending_points_ID
    rownames(find_path_distances) <- starting_points_ID

    # Matches starting point with ending point
    greedy_solution_translated <- greedy_match(find_path_distances)

    # We get the edges between the starting points and the ending points
    keep_edges     <- vector(mode = "list", length = npoints)
    for (j in seq_len_npoints){
      keep_edges[[j]] <- as.numeric(
        shortest_paths(
          graph  = net_result_loop,
          from   = greedy_solution_translated[j,1],
          to     = greedy_solution_translated[j,2],
          output = "epath",
          mode   = "out"
        )$epath[[1]]
      )
    }

    table_keep_edges     <- table(unlist(keep_edges))
    keep_edges_dataframe <- data.frame(
      link      = as.numeric(names(table_keep_edges)),
      frequency = as.numeric(table_keep_edges)
    )

    merged <- merge(x     = initial_w_matrix,
                    y     = keep_edges_dataframe,
                    by    = "link",
                    all.x = TRUE)
    merged$frequency[is.na(merged$frequency)] <- 0
    merged$frequency                          <- merged$frequency / population

    if (i == 1){
      utilities           <- data.frame("link"  = merged$link)
      utilities$edge_cost <- merged$frequency + merged$weight
    } else {
      utilities$edge_cost <- utilities$edge_cost + merged$frequency
    }

    net_result_loop     <- set_edge_attr(graph = net_result_loop,
                                         name  = "weight",
                                         index = E(net_result_loop),
                                         value = utilities$edge_cost)

    congestion                    <- utilities
    names(congestion)[2]          <- "c_l"
    congestion$u_l                <- congestion$c_l - merged$weight
    congestion$u_l_population     <- congestion$u_l * population
    congestion$u_l_perc           <- congestion$u_l_population / (npoints * i)
    congestion$u_l_perc_times_c_l <- congestion$u_l_perc * congestion$c_l

    list_of_congestion[[i]] <- congestion
    list_of_edges_freq[[i]] <- keep_edges_dataframe

    if (all_networks){

      E_net_result <- E(net_result_loop)

      net_result_congestion <- set_edge_attr(graph = net_result_loop,
                                             name  = "c_l",
                                             index = E_net_result,
                                             value = congestion$c_l
      )

      net_result_congestion <- set_edge_attr(graph = net_result_congestion,
                                             name  = "u_l",
                                             index = E_net_result,
                                             value = congestion$u_l
      )

      net_result_congestion <- set_edge_attr(graph = net_result_congestion,
                                             name  = "u_l_population",
                                             index = E_net_result,
                                             value = congestion$u_l_population
      )

      net_result_congestion <- set_edge_attr(graph = net_result_congestion,
                                             name  = "u_l_perc",
                                             index = E_net_result,
                                             value = congestion$u_l_perc
      )

      net_result_congestion <- set_edge_attr(graph = net_result_congestion,
                                             name  = "u_l_perc_times_c_l",
                                             index = E_net_result,
                                             value =
                                               congestion$u_l_perc_times_c_l
      )

      list_of_networks[[i]] <- net_result_congestion
    }

  }

  E_net_result <- E(net_result)

  net_result_congestion <- set_edge_attr(graph = net_result,
                                         name  = "c_l",
                                         index = E_net_result,
                                         value = congestion$c_l
  )

  net_result_congestion <- set_edge_attr(graph = net_result_congestion,
                                         name  = "u_l",
                                         index = E_net_result,
                                         value = congestion$u_l
  )

  congestion$u_l_perc_times_c_l <- congestion$u_l_perc * congestion$c_l

  net_result_congestion <- set_edge_attr(graph = net_result_congestion,
                                         name  = "u_l_population",
                                         index = E_net_result,
                                         value = congestion$u_l_population
  )

  net_result_congestion <- set_edge_attr(graph = net_result_congestion,
                                         name  = "u_l_perc",
                                         index = E_net_result,
                                         value = congestion$u_l_perc
  )

  net_result_congestion <- set_edge_attr(graph = net_result_congestion,
                                         name  = "u_l_perc_times_c_l",
                                         index = E_net_result,
                                         value = congestion$u_l_perc_times_c_l
  )

  solution_edges <- st_as_sf(net_result_congestion, "edges")

  metrics_df <- data.frame("niters"                    = seq_int_niters,
                           "rmse_u_l_perc"             = NA,
                           "rmse_u_l_perc_times_c_l"   = NA,
                           "rmse_c_l"                  = NA,
                           "sum_c_l"                   = NA
  )

  metrics_df$sum_c_l[1] <- sum(list_of_congestion[[1]]$c_l)

  metrics_df[,1] <- seq_int_niters
  for (i in seq_int_niters[-1]){
    metrics_df$rmse_u_l_perc[i]           <-
      sqrt(mean(
        (list_of_congestion[[i-1]]$u_l_perc -
           list_of_congestion[[i]]$u_l_perc)^2
      ))

    metrics_df$rmse_u_l_perc_times_c_l[i] <-
      sqrt(mean(
        (list_of_congestion[[i-1]]$u_l_perc_times_c_l -
           list_of_congestion[[i]]$u_l_perc_times_c_l)^2
      ))

    metrics_df$rmse_c_l[i] <-
      sqrt(mean(
        (list_of_congestion[[i-1]]$c_l -
           list_of_congestion[[i]]$c_l)^2
      ))

    metrics_df$sum_c_l[i] <- sum(list_of_congestion[[i]]$c_l)
  }

  result <- list(
    solution_edges        = solution_edges,
    net_result_congestion = net_result_congestion,
    list_of_edges_freq    = list_of_edges_freq,
    list_of_congestion    = list_of_congestion,
    metrics_df            = metrics_df,
    origin_points         = random_points_1,
    destination_points    = random_points_2,
    list_of_networks      = list_of_networks
  )

  class(result) <- c("EGCorrs", class(result))

  return(result)

}
