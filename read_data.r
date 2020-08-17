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

