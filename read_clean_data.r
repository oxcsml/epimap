Cinfer_mean <- read.csv('data/Cinfer_mean.csv')
colnames(Cinfer_mean) <- dates[1:dim(Cinfer_sample)[2]]
rownames(Cinfer_mean) <- areas
Nsample = 10
Cinfer_sample = list()
Clean_sample = list()
for (i in 1:Nsample) {
  Cinfer_sample[[i]] <- read.csv(paste('data/Cinfer_sample',i,'.csv',sep=''))
  colnames(Cinfer_sample[[i]]) <- dates[1:dim(Cinfer_sample[[i]])[2]]
  rownames(Cinfer_sample[[i]]) <- areas
  Clean_sample[[i]] <- read.csv('data/Clean_sample.csv')
  colnames(Clean_sample[[i]]) <- dates[1:dim(Clean_sample[[i]])[2]]
  rownames(Clean_sample[[i]]) <- areas
}

