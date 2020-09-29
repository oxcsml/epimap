library(geosphere)

# Case time series in GB areas
uk_cases <- read.csv("data/uk_cases.csv")
ind <- sapply(uk_cases[,2], function(s)
    !(s %in% c('Outside Wales','Unknown','...17','...18'))
)
uk_cases <- uk_cases[ind,]

N <- nrow(uk_cases) # 348 LTLA in England, NHS health boards Scotland, X in Wales

areas <- uk_cases[,2]
quoted_areas <- sapply(areas,function(s) paste('"',s,'"',sep=''))

###############################################################################

# Assumes that metadata contains all areas in uk_cases and area names match perfectly
metadata <- read.csv("data/metadata.csv")
ind <- sapply(metadata[,1], function(s)
    s %in% areas
)
metadata <- metadata[ind,]

geoloc <- matrix(0, N, 2)
geodist <- matrix(0, N, N)
population <- rep(0.0, N)

# match up areas in uk_cases and metadata
meta_areas <- metadata$AREA
longitudes <- metadata$LONG
latitudes <- metadata$LAT

for (i in 1:N) {
  area <- areas[i]
  j <- grep(sprintf('^%s$',area), meta_areas)
  l <- length(j)
  if (l >= 1) {                 # just use first match!!
    if (l > 1) {
      print(sprintf("Area: %s",area))
      print(sprintf("Matched meta_areas: %s", paste(meta_areas[j],collapse=', ')))
      print(sprintf("Using: %s, %s",meta_areas[j[l]],metadata$CODE[j[l]]))
    }
    geoloc[i, 1] = longitudes[j[l]]
    geoloc[i, 2] = latitudes[j[l]]
    population[i] = metadata$POPULATION[j[l]]
  } else {
    print(sprintf("Cannot find area '%s'",area))
    for (r in 1:length(meta_areas)) {
      if (length(grep(meta_areas[r], area))>0) {
        geoloc[i, 1] = longitudes[r]
        geoloc[i, 2] = latitudes[r]
        population[i] = metadata$POPULATION[r]
        print(sprintf("...found area '%s'",meta_areas[r]))
      }
    }
  }
}

write.csv(
  data.frame(
    area=quoted_areas, 
    longitude=geoloc[,1], 
    latitude=geoloc[,2], 
    population=population
  ),'data/areas.csv',row.names=FALSE,quote=FALSE)



#########################################################################

write.matrix.csv <- function(m,filename) {
  colnames(m) <- quoted_areas
  rownames(m) <- quoted_areas
  write.csv(m,filename,quote=FALSE)
}
  
# compute distances between areas. Straight line, not actual travel distance.
for (i in 1:N) {
  for (j in i:N) {
    # distance between two points on an ellipsoid (default is WGS84 ellipsoid), in units of 100km
    geodist[i, j] = distGeo(geoloc[i, 1:2], geoloc[j, 1:2]) / 100000
    geodist[j, i] = geodist[i, j]
  }
}
write.matrix.csv(geodist,'data/distances.csv')

