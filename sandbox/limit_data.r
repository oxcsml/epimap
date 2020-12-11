limit_data_by_distance <- function(env,area,distance) { 
  ind <- env$geodist[area,] < distance
  env$N <- sum(ind)
  env$areas <- env$areas[ind]
  env$quoted_areas <- env$quoted_areas[ind]
  
  env$sparse_region <- env$sparse_region[ind,]
  region_ind <- colSums(env$sparse_region != 0) !=0
  env$sparse_region <- env$sparse_region[,region_ind]
  env$N_region <- ncol(env$sparse_region)
  if (is.null(env$N_region)) {
    env$sparse_region = matrix(env$sparse_region)
    env$N_region = 1
  }
  env$quoted_regions <- env$quoted_regions[region_ind]

  env$geoloc <- env$geoloc[ind,]
  env$geodist <- env$geodist[ind,ind]
  env$AllCount <- env$AllCount[ind,]
  env$Clean_latent <- env$Clean_latent[ind,]
  env$Clean_recon <- env$Clean_recon[ind,]
  env$radiation_flux <- env$radiation_flux[ind,ind,]
  env$traffic_flux <- env$traffic_flux[ind,ind,]
  for (i in 1:env$N) {
    for (f in 1:dim(env$radiation_flux)[3]) {
      x <- env$radiation_flux[i,,f]
      env$radiation_flux[i,,f] <- x / sum(x)
    }
    for (f in 1:dim(env$traffic_flux)[3]) {
      x <- env$traffic_flux[i,,f]
      env$traffic_flux[i,,f] <- x / sum(x)
    }
  }
}

limit_data_multi <- function(env,areas_distances) { 
  ind <- apply(sapply(areas_distances,
    function(ad) env$geodist[ad$area,] < ad$distance
  ),c(1),any)
  env$N <- sum(ind)
  env$areas <- env$areas[ind]
  env$quoted_areas <- env$quoted_areas[ind]
  
  env$sparse_region <- env$sparse_region[ind,]
  region_ind <- colSums(env$sparse_region != 0) !=0
  env$sparse_region <- env$sparse_region[,region_ind]
  env$N_region <- ncol(env$sparse_region)
  if (is.null(env$N_region)) {
    env$sparse_region = matrix(env$sparse_region)
    env$N_region = 1
  }
  env$quoted_regions <- env$quoted_regions[region_ind]

  env$geoloc <- env$geoloc[ind,]
  env$geodist <- env$geodist[ind,ind]
  env$AllCount <- env$AllCount[ind,]
  env$Clean_latent <- env$Clean_latent[ind,]
  env$Clean_recon <- env$Clean_recon[ind,]
  env$radiation_flux <- env$radiation_flux[ind,ind,]
  env$traffic_flux <- env$traffic_flux[ind,ind,]
  for (i in 1:env$N) {
    for (f in 1:dim(env$radiation_flux)[3]) {
      x <- env$radiation_flux[i,,f]
      env$radiation_flux[i,,f] <- x / sum(x)
    }
    for (f in 1:dim(env$traffic_flux)[3]) {
      x <- env$traffic_flux[i,,f]
      env$traffic_flux[i,,f] <- x / sum(x)
    }
  }
}
