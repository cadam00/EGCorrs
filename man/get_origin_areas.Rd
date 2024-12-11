\name{get_origin_areas}

\alias{get_origin_areas}

\title{Example mask shapefile}

\description{
Example mask shapefile.
}

\usage{
get_origin_areas()
}

\details{
Example of input \code{origin_areas} used for functions. It is a sample
from Natura 2000 areas.
}

\value{
\code{sfc_MULTIPOLYGON} and \code{sfc} object.
}

\references{
European Environment Agency. (2023). Natura 2000 data - the European network of
protected sites. European Environment Agency. Retrieved from
\href{https://www.eea.europa.eu/en}{https://www.eea.europa.eu/en}
Last Access: 16/10/2024.
}

\examples{
origin_areas <- get_origin_areas()
plot(origin_areas)
}
