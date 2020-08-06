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

Count <- uk_cases[,3:ncol(uk_cases)]
dates <- as.Date(colnames(uk_cases)[3:ncol(uk_cases)], format='X%Y.%m.%d')

rownames(Count) <- quoted_areas
colnames(Count) <- dates

write.csv(Count,'data/counts.csv',quote=FALSE)

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

#########################################################################

# compute fluxes for radiation model.
# relax radiation model by adding noise to distance with sd 10km, and bootstrapping

numrep <- 100000

for (radiation_length_scale in c(.1,.2,.5)) {
  flux <- array(0, dim=c(N, N))
  flux_var <-  array(0, dim=c(N, N))
  for (r in 1:numrep) {
    ff <- matrix(0, N, N)
    for (i in 1:N) {
      mi = population[i]
      if (i==1) {
        jj <- 2:N
      } else if (i==N) {
        jj <- 1:(N-1)
      } else {
        jj = c(1:(i-1),(i+1):N)
      }
      js = jj[order(geodist[i,jj]+rnorm(N-1,sd=radiation_length_scale))]
      nj = population[js]
      sij = cumsum(c(0.0,nj))[1:(N-1)]
      ff[i,js] = mi*nj/((mi+sij)*(mi+nj+sij))
      ff[i,i] = 1.0 - sum(ff[i,js]) # this statement weird, as fluxes don't sum to 1, unlike as claimed in Simini paper
    }
    flux = flux + ff
    flux_var = flux_var + ff^2
  }
  flux = flux / numrep
  flux_sd = (flux_var - numrep*(flux^2))/(numrep-1)
  flux_sd[flux_sd<0]=0
  flux_sd = sqrt(flux_sd)/sqrt(numrep)

  write.matrix.csv(flux,sprintf('data/radiation_flux_ls=%1.1f.csv',radiation_length_scale))
  write.matrix.csv(flux_sd,sprintf('data/radiation_flux_sd_ls=%1.1f.csv',radiation_length_scale))
}
