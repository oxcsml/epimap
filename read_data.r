infprofile <- read.csv("data/serial_interval.csv")$fit
D <- length(infprofile)

df <- read.csv("data/areas.csv",row.names=1)
N <- nrow(df)
areas <- rownames(df)
quoted_areas <- sapply(areas,function(s) paste('"',s,'"',sep=''))
geoloc <- df[,1:2]
population <- df[,3]

Count <- read.csv("data/counts.csv",row.names=1)
Tall <- ncol(Count)
dates <- as.Date(colnames(Count), format='X%Y.%m.%d')
colnames(Count) <- dates

geodist <- read.csv("data/distances.csv",row.names=1)
colnames(geodist) <- areas

radiation_length_scales <- c(.1,.2,.5)
flux <- array(0,dim=c(N,N,length(radiation_length_scales)))
for (i in 1:length(radiation_length_scales)) {
  ls <- radiation_length_scales[i]
  df <- data.matrix(read.csv(sprintf('data/radiation_flux_ls=%1.1f.csv',ls),row.names=1))
  flux[,,i] <- df
}
colnames(flux) <- areas
rownames(flux) <- areas
dimnames(flux)[[3]] <- radiation_length_scales

