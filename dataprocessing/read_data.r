##########################################################################
##########################################################################
Rmap_read_data = function(env) { 
  with(env,{

    readdata = function(filename,...) {
      read.csv(sprintf('%s/%s.csv',opt$data_directory,filename),...)
    }
    readclean = function(filename,...) {
      read.csv(sprintf('%s/%s.csv',opt$clean_directory,filename),...)
    }

    #########################################################
    infprofile <- readdata("serial_interval")$fit
    Tip <- 30
    infprofile <- infprofile[1:Tip]
    infprofile <- infprofile/sum(infprofile)
    D <- length(infprofile)

    Tdp <- 14
    Adp <- 5.0
    Bdp <- 1.0
    testdelayprofile <- pgamma(1:Tdp,shape=Adp,rate=Bdp)
    testdelayprofile <- testdelayprofile/testdelayprofile[Tdp]
    testdelayprofile <- testdelayprofile - c(0.0,testdelayprofile[1:(Tdp-1)])


    df <- readdata("areas",row.names=1)
    N <- nrow(df)
    areas <- rownames(df)
    quoted_areas <- sapply(areas,function(s) paste('"',s,'"',sep=''))
    geoloc <- df[,1:2]
    population <- df[,3]

    geodist <- readdata("distances",row.names=1)
    colnames(geodist) <- areas

    # Use counts from uk_cases in case updated
    uk_cases <- readdata("uk_cases")
    ind <- sapply(uk_cases[,2], function(s)
        !(s %in% c('Outside Wales','Unknown','...17','...18'))
    )
    uk_cases <- uk_cases[ind,]
    AllCount <- uk_cases[,3:ncol(uk_cases)]
    Tall <- ncol(AllCount)
    dates <- as.Date(colnames(AllCount), format='X%Y.%m.%d')
    colnames(AllCount) <- dates
    rownames(AllCount) <- areas

    #########################################################
    radiation_length_scales <- c(.1,.2,.5)
    radiation_flux <- array(0,dim=c(N,N,length(radiation_length_scales)))
    for (i in 1:length(radiation_length_scales)) {
      ls <- radiation_length_scales[i]
      df <- data.matrix(readdata(sprintf('radiation_flux_ls=%1.1f',ls),row.names=1))
      radiation_flux[,,i] <- df
    }
    colnames(radiation_flux) <- areas
    rownames(radiation_flux) <- areas
    dimnames(radiation_flux)[[3]] <- radiation_length_scales

    traffic_flux <- array(0, dim=c(N,N,2))
    df <- data.matrix(readdata('traffic_flux_row-normed', row.names=1))
    traffic_flux[,,1] <- df
    df <- data.matrix(readdata('traffic_flux_transpose_row-normed', row.names=1))
    traffic_flux[,,2] <- df
    colnames(traffic_flux) <- areas
    rownames(traffic_flux) <- areas

    #########################################################
    # TODO: Second block is redundant here?
    if (opt$observation_data == 'cleaned_latent_sample' ||
        opt$observation_data == 'cleaned_recon_sample') {
      sample_id = opt$cleaned_sample_id
      Clean_latent <- readclean(paste('Clatent_sample',sample_id,sep=''), row.names=1)
      Clean_recon <- readclean(paste('Crecon_sample',sample_id,sep=''), row.names=1)
      print(paste('Using samples from Clatent_sample',sample_id,'.csv',sep=''))
    } else {
      sample_id = 'mean'
      Clean_latent <- readclean('Clatent_mean', row.names=1)
      Clean_recon <- readclean('Crecon_median', row.names=1)
      # placeholder if not using cleaned data
    }

  })
  env
}# Rmap_read_data
##########################################################################
##########################################################################
