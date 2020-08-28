
radiation_length_scales <- c(.1,.2,.5)
radiation_flux <- array(0,dim=c(N,N,length(radiation_length_scales)))
for (i in 1:length(radiation_length_scales)) {
  ls <- radiation_length_scales[i]
  df <- data.matrix(read.csv(sprintf('data/radiation_flux_ls=%1.1f.csv',ls),row.names=1))
  radiation_flux[,,i] <- df
}
colnames(radiation_flux) <- areas
rownames(radiation_flux) <- areas
dimnames(radiation_flux)[[3]] <- radiation_length_scales

traffic_flux <- array(0, dim=c(N,N,1))
df <- data.matrix(read.csv('data/traffic_flux.csv', row.names=1))
traffic_flux[,,1] <- df
colnames(traffic_flux) <- areas
rownames(traffic_flux) <- areas
