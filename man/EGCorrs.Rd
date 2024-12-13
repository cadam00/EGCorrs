\name{EGCorrs}
\alias{EGCorrs}
\title{Species Corridors Estimation}
\description{
Sea Currents To Evolutionary Game Theoretic Model for Species Corridors
Estimation
}
\usage{
EGCorrs(component_u, component_v, origin_areas, destination_areas,
         npoints = 40L, lambda = 1, niters = 100L, k_neighbors = 7L,
         nearest_grid_nodes = 4L, mask_shapefile = NULL, all_networks = FALSE,
         progBar = TRUE)
}

\arguments{
\item{component_u}{
\code{SpatRaster} object with horizontal direction.
}
\item{component_v}{
\code{SpatRaster} object with vertical direction.
}
\item{origin_areas}{
\code{sf} object with origin areas, as \code{POLYGON} or \code{MULTIPOLYGON}
rows.
}
\item{destination_areas}{
\code{sf} object with destination areas, as \code{POLYGON} or
\code{MULTIPOLYGON} rows.
}
\item{npoints}{
\code{integer} object with the number of random points for origin and
destination areas. It must be positive. The default is \code{40L}.
}
\item{lambda}{
\code{numeric} object with the aggregation parameter \eqn{\lambda}. The default
value is \code{1}.
}
\item{niters}{
\code{integer} object with the number of points for origin and destination
areas. It must be positive. The default value is \code{100L}.
}
\item{k_neighbors}{
\code{integer} object indicating the number of \eqn{k} of the nearest neighbors
to use from the components of horizontal and vertical directions. The default
value is \code{7L} nearest neighbors.
}
\item{nearest_grid_nodes}{
\code{integer} object indicating the number of \eqn{k} of the nearest neighbors
to use from the random points inside origin and destination areas. The default
value is \code{4L} nearest neighbors.
}
\item{mask_shapefile}{
Optional \code{sf} or SpatVector object for masking \code{component_u} and
\code{component_v} \code{SpatRaster} objects. The default value is \code{NULL},
indicating that no mask is performed.
}
\item{all_networks}{
\code{logical} object indicating the return of \code{sfnetwork} per iteration.
}
\item{progBar}{
\code{logical} object indicating whether to print progress bar or not. The
default is printing (\code{progBar=TRUE}).
}
}

\details{
Sea currents data of the input are transformed into weighted directed graph
connectivity data. Each centroid of \code{component_u} and \code{component_v} is
considered as a separate graph node. Except from these points, random
\code{npoints} random origin points are generated inside the \code{origin_areas}
polygons/multipolygons and \code{npoints} random destination points are
generated inside the \code{destination_areas} polygons/multipolygons. These
points are equally or near equally distributed between each group of points, in
a way that totally \code{npoints} random points are contained inside
\code{origin_areas} and \code{npoints} random points are contained inside
\code{destination_areas}. The directions produced from these points are
calculated, multiplied by \eqn{\lambda} and finally minmax-scaled in
\eqn{[0-1]}. Thereafter, shortest paths between origin and destination points
are computed \code{niters} times, following a evolutionary game theoretic model.

Extent, resolution and coordinate system of both must have the same for both
\code{component_u} and \code{component_v}; otherwise an error stops the
function's execution. In case of mismatches between the area of these two
components and \code{origin_areas} or \code{destination_areas}, then the
intersection of these areas and the components is used instead of these areas.
In case that a \code{mask_shapefile} is provided, then it is internally assured
that it has the same resolution or coordinates with \code{component_u} and a
warning is prompted. Reproducibility of the results could be guranteed by the
use of \code{base::set.seed} before using the function.
}

\value{
\code{SeaGraphs} object, which is also a \code{list} object containing the
following elements:
\itemize{
\item{solution_edges}: \code{sf} object representing the final edges of the
solution after \code{niters} iterations.
\item{net_result_congestion}: \code{sfnetwork} object representing the final
solution after \code{niters} iterations. Each row has the following items for
each graph edge/link \eqn{l}:
\itemize{
\item{from}: edge start for a graph vertex ID.
\item{to}: edge end for a graph vertex ID.
\item{c_l}: value of cost function \eqn{c_l}.
\item{u_l}: value of utility \eqn{u_l}.
\item{u_l_population}: value of utility compared to the total population of the
game (where \eqn{population = niters * npoints}).
\item{u_l_perc}: value of u_l_population normalized by the population of the
current iteration (\eqn{u_{l_{_{perc}}}}).
\item{u_l_perc_times_c_l}: u_l_perc * c_l.
}
\item{list_of_edges_freq}: \code{list} object containing \code{data.frame}
objects with the the frequency that shortest path corridors pass through each
graph vertex, for each iteration.
\item{list_of_congestion}: \code{list} object containing \code{data.frame}
objects with congenstion metrics per graph edge for each one of the
\code{niters} iterations. Each row has a number indicating the graph edge/link
ID, according to the column link. All the other columns are as in
net_result_congestion output.
\item{metrics_df}: \code{data.frame} object with summarization of
list_of_congestion values accross, where:
\itemize{
\item{niters}: number of the corresponding iteration.
\item{rmse_u_l_perc}: \eqn{RMSE(u_{l_{_{perc}}})}.
\item{rmse_u_l_perc_times_c_l}: \eqn{RMSE(u_{l_{_{perc}}} * c_l)}.
\item{rmse_c_l}: \eqn{RMSE(c_l)}.
\item{sum_c_l}: \eqn{\sum _l{c_l}}.
}
\item{origin_points}: \code{sf} \code{POINT} object representing the origin
points.
\item{destination_points}: \code{sf} \code{POINT} object representing the
destination points.
\item{list_of_networks}: If \code{all_networks = TRUE}, then it is a \code{list}
containing net_result_congestion item for each iteration. Else, it is
\code{NULL}.
}
}

\examples{
set.seed(42)
component_u       <- get_component_u()
component_v       <- get_component_v()
origin_areas      <- get_origin_areas()
destination_areas <- get_destination_areas()
# Warnings are thrown because of reassuring that components u/v and
# origin/destination areas are intersecting
out               <- EGCorrs(component_u, component_v, origin_areas,
                         destination_areas, npoints = 1,  niters = 10,
                         progBar=FALSE)
## Warning messages:
## 1: attribute variables are assumed to be spatially constant throughout all
## geometries
## 2: attribute variables are assumed to be spatially constant throughout all
## geometries
}
