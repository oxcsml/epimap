source("dataprocessing/read_data.r")
source("mapping/epimap.r")
source("dataprocessing/process_radiation_fluxes.r")

opt = Rmap_options()

readdata = function(filename,...) {
    read.csv(sprintf('%s/%s.csv',opt$data_directory,filename),...)
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

  write.matrix.csv(val$flux,sprintf('data/radiation_flux_ls=%1.1f_2.csv',radiation_length_scale))
  write.matrix.csv(val$flux_sd,sprintf('data/radiation_flux_sd_ls=%1.1f_2.csv',radiation_length_scale))
}
