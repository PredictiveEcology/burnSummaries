## One polygon inside each of the cells `cells` of raster `r`, with fire years `years`, shaped as
## fireregimetools' loaders return them. Each covers its cell's centre and no other cell's.
cellPolys <- function(r, cells, years) {
  xy <- terra::xyFromCell(r, cells)
  h <- 0.4 * terra::res(r)
  wkt <- sprintf(
    "POLYGON ((%1$f %3$f, %2$f %3$f, %2$f %4$f, %1$f %4$f, %1$f %3$f))",
    xy[, 1] - h[1], xy[, 1] + h[1], xy[, 2] - h[2], xy[, 2] + h[2]
  )
  v <- terra::vect(wkt, crs = terra::crs(r))
  v$YEAR <- as.integer(years)
  v$SIZE_HA <- rep(1, length(years))
  v
}
