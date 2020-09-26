limit_data_by_distance <- function(env,area,distance) { 
  ind <- env$geodist[area,] < distance
  env$N <- sum(ind)
  env$areas <- env$areas[ind]
  env$quoted_areas <- env$quoted_areas[ind]
  env$geoloc <- env$geoloc[ind,]
  env$geodist <- env$geodist[ind,ind]
  env$Count <- env$Count[ind,]
  env$AllCount <- env$AllCount[ind,]
  env$Clean_latent <- env$Clean_latent[ind,]
  env$Clean_recon <- env$Clean_recon[ind,]
  env$radiation_flux <- env$radiation_flux[ind,ind,]
  for (i in 1:env$N) {
    for (f in 1:dim(env$radiation_flux)[3]) {
      x <- env$radiation_flux[i,,f]
      env$radiation_flux[i,,f] <- x / sum(x)
    }
  }
}
  
