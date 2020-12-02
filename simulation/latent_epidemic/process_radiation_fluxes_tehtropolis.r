source("mapping/epimap.r")
source("dataprocessing/process_radiation_fluxes.r")

readdata = function(filename,...) {
    read.csv(sprintf('%s/%s.csv',"simulation/latent_epidemic/data",filename),...)
}

df <- readdata("areas",row.names=1)
N <- nrow(df)
areas <- rownames(df)
quoted_areas <- sapply(areas,function(s) paste('"',s,'"',sep=''))
geoloc <- df[,1:2]
population <- df[,3]
population <- as.double(population)

geodist <- readdata("distances",row.names=1)
colnames(geodist) <- areas

for (radiation_length_scale in c(.1,.2,.5)) {
  numrep <- 100000
  val = compute_radiation_fluxes(numrep, radiation_length_scale, population, geodist)

  write.csv(val$flux,sprintf('simulation/latent_epidemic/data/radiation_flux_ls=%1.1f.csv',radiation_length_scale))
  write.csv(val$flux_sd,sprintf('simulation/latent_epidemic/data/radiation_flux_sd_ls=%1.1f.csv',radiation_length_scale))
}
