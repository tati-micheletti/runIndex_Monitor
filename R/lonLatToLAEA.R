#' lon/lat (degrees) -> ETRS89-LAEA Europe (EPSG:3035) easting/northing (m), closed form (Snyder 1987; EPSG Guidance Note 7-2)
#'
#' A copy of dataPrep_Monitor's `lonLatToLAEA()` (with its GRS80 constants inside) so this module does not depend on another one; used by studyArea.R. Pure R, NO GDAL/PROJ:
#' the EPSG:3035 axis order of terra/sf transformations has come out swapped on EVE (DECISIONS.md 2026-10-05, 2026-10-08).
#' @param lon,lat Numeric vectors (degrees). @return Matrix with columns x (easting), y (northing).
lonLatToLAEA <- function(lon, lat) {
  a <- 6378137; f <- 1 / 298.257222101; e2 <- 2 * f - f^2; e <- sqrt(e2)
  qf <- function(phi) { s <- sin(phi); (1 - e2) * (s / (1 - e2 * s^2) - (1 / (2 * e)) * log((1 - e * s) / (1 + e * s))) }
  phi <- lat * pi / 180; lam <- lon * pi / 180; phi0 <- 52 * pi / 180; lam0 <- 10 * pi / 180
  q <- qf(phi); qp <- qf(pi / 2); q0 <- qf(phi0)
  Rq <- a * sqrt(qp / 2); beta <- asin(q / qp); beta0 <- asin(q0 / qp)
  m1 <- cos(phi0) / sqrt(1 - e2 * sin(phi0)^2); D <- a * m1 / (Rq * cos(beta0))
  B <- Rq * sqrt(2 / (1 + sin(beta0) * sin(beta) + cos(beta0) * cos(beta) * cos(lam - lam0)))
  cbind(x = 4321000 + B * D * cos(beta) * sin(lam - lam0),
        y = 3210000 + (B / D) * (cos(beta0) * sin(beta) - sin(beta0) * cos(beta) * cos(lam - lam0)))
}
