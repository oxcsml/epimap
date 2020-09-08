limit_data_by_distance <- function(area,distance) {
  ind <- geodist[area,] < distance
  N <<- sum(ind)
  areas <<- areas[ind]
  quoted_areas <<- quoted_areas[ind]
  geoloc <<- geoloc[ind,]
  geodist <<- geodist[ind,ind]
  Count <<- Count[ind,]
  flux <- radiation_flux[ind,ind,]
  for (i in 1:N) {
    for (f in 1:dim(flux)[3]) {
      x <- flux[i,,f]
      flux[i,,f] <- x / sum(x)
    }
  }
  radiation_flux <<- flux
} 
  
