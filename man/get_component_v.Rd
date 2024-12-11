\name{get_component_v}

\alias{get_component_v}

\title{Example component v}

\description{
Example horizontal direction raster.
}

\usage{
get_component_v()
}

\details{
Example of input \code{component_v} used for functions. It is a cropped area of
the Aegean Sea.
}

\value{
\code{SpatRaster} object.
}

\references{
Clementi, E., Aydogdu, A., Goglio, A. C., Pistoia, J., Escudier, R., Drudi, M.,
Grandi, A., Mariani, A., Lyubartsev, V., Lecci, R., Cretí, S., Coppini, G.,
Masina, S., & Pinardi, N. (2021). Mediterranean Sea Physical Analysis and
Forecast (CMEMS MED-Currents, EAS6 system) (Version 1) [Data set]. Copernicus
Monitoring Environment Marine Service (CMEMS).
\doi{10.25423/CMCC/MEDSEA_ANALYSISFORECAST_PHY_006_013_EAS8}.
Last Access: 16/10/2024.
}

\seealso{
\code{
\link[EGCorrs]{get_component_u}
}
}

\examples{
component_v <- get_component_v()
terra::plot(component_v)
}
