\name{latlong2km}
\alias{latlong2km}
\title{Convert latitude and longitude to Northings and Eastings}
\usage{
  latlong2km(lon, lat, lon0 = sum(range(lon))/2,
    lat0 = sum(range(lat))/2)
}
\arguments{
  \item{lon}{longitude}

  \item{lat}{latitude}

  \item{lon0}{longitude reference point (defaults to mean
  longitude)}

  \item{lat0}{latitude reference point (defaults to mean
  latitude)}
}
\value{
  list with elements \code{km.e} and \code{km.n}.
}
\description{
  Convert longitude and latitude co-ordinates to kilometres
  west-east and south-north from axes through
  (\code{lon0},\code{lat0}) using the "spherical law of
  cosines".
}
\author{
  Simon N. Wood
}

