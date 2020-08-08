infprofile <- read.csv("data/serial_interval.csv")$fit
D <- length(infprofile)

df <- read.csv("data/areas.csv",row.names=1)
N <- nrow(df)
areas <- rownames(df)
quoted_areas <- sapply(areas,function(s) paste('"',s,'"',sep=''))
geoloc <- df[,1:2]
population <- df[,3]

geodist <- read.csv("data/distances.csv",row.names=1)
colnames(geodist) <- areas

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

# Use counts from uk_cases in case updated

uk_cases <- read.csv("data/uk_cases.csv")
ind <- sapply(uk_cases[,2], function(s)
    !(s %in% c('Outside Wales','Unknown','...17','...18'))
)
uk_cases <- uk_cases[ind,]
Count <- uk_cases[,3:ncol(uk_cases)]
Tall <- ncol(Count)
dates <- as.Date(colnames(Count), format='X%Y.%m.%d')
colnames(Count) <- dates
rownames(Count) <- areas

