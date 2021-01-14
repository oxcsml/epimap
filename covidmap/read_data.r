library("rjson")

covidmap_read_data = function(env) { 
  with(env,{

    readdata = function(filename,...) {
      read.csv(sprintf('%s/%s.csv',opt$data_directory,filename),...)
    }
    readsinglearea = function(filename,...) {
      read.csv(sprintf('%s/%s.csv',opt$singlearea_directory,filename),...)
    }

    #########################################################
    infprofile <- readdata("serial_interval")$fit
    Tip <- 30
    infprofile <- infprofile[1:Tip]
    infprofile <- infprofile/sum(infprofile)
    D <- length(infprofile)

    Tdp <- 21
    presymptomdays <- 2
    Tdpnz <- Tdp - presymptomdays
    Adp <- 5.8
    Bdp <- 0.948
    testdelayprofile <- pgamma(1:Tdpnz,shape=Adp,rate=Bdp)
    testdelayprofile <- testdelayprofile/testdelayprofile[Tdpnz]
    testdelayprofile <- testdelayprofile - c(0.0,testdelayprofile[1:(Tdpnz-1)])
    testdelayprofile <- c(rep(0, presymptomdays), testdelayprofile)

    # data about areas
    df <- readdata("areas",row.names=1)
    N <- nrow(df)
    areas <- rownames(df)
    quoted_areas <- sapply(areas,function(s) paste('"',s,'"',sep=''))
    geoloc <- df[,1:2]
    population <- df[,3]

    # NHS region data
    sparse_region <- df[,4:12]
    N_region <- ncol(sparse_region)
    region_df <- readdata("nhs_regions", row.names=1)
    regions <- rownames(region_df)
    quoted_regions <- sapply(regions,function(s) paste('"',s,'"',sep=''))

    geodist <- readdata("distances",row.names=1)
    colnames(geodist) <- areas

    # read data about how to split UK into regions to model for full model
    js = fromJSON(file="data/region-groupings.json") ### TODO FIX
    inferred_region = sparse_region
    modelled_region = 0*sparse_region
    for (i in 1:9) {
      modelled_region[js[[i]],i]=1
    }
    stopifnot(all(inferred_region<=modelled_region))

    # Use counts from uk_cases in case updated
    cases <- readdata("cases")
    ind <- sapply(cases[,2], function(s)
        !(s %in% c('Outside Wales','Unknown','...17','...18'))
    )
    cases <- cases[ind,]
    AllCount <- cases[,3:ncol(cases)]
    Tall <- ncol(AllCount)
    dates <- as.Date(colnames(AllCount), format='X%Y.%m.%d')
    colnames(AllCount) <- dates
    rownames(AllCount) <- areas

    if (!identical(opt$approximation, "singlearea")) {
      if (opt$singlearea_sample_id>0) {
        sample_id = opt$singlearea_sample_id
        Clean_latent = readsinglearea(paste('Clatent_sample',sample_id,sep=''), 
          row.names=1)
        Clean_recon = readsinglearea(paste('Crecon_sample',sample_id,sep=''), 
          row.names=1)
        print(paste(
         'Using samples from ',
         opt$singlearea_directory,'/Clatent_sample',sample_id,'.csv',
         sep=''
        ))
      } else {
        Clean_latent <- readsinglearea('Clatent_median', row.names=1)
        Clean_recon <- readsinglearea('Crecon_median', row.names=1)
        print(paste(
          'Using samples from ',
          opt$singlearea_directory,'/Clatent_median',
          sep=''
        ))
      }
    }

    #########################################################
    # Radiation fluxes have also been tried, but found to not be significantly different from traffic models.
    # radiation_length_scales <- c(.1,.2,.5)
    # radiation_flux <- array(0,dim=c(N,N,length(radiation_length_scales)))
    # for (i in 1:length(radiation_length_scales)) {
    #   ls <- radiation_length_scales[i]
    #   df <- data.matrix(readdata(sprintf('radiation_flux_ls=%1.1f',ls),row.names=1))
    #   radiation_flux[,,i] <- df
    # }
    # colnames(radiation_flux) <- areas
    # rownames(radiation_flux) <- areas
    # dimnames(radiation_flux)[[3]] <- radiation_length_scales

    traffic_flux <- array(0, dim=c(N,N,2))
    traffic_flux[,,1] <- data.matrix(readdata('traffic_flux_row-normed', row.names=1))
    traffic_flux[,,2] <- data.matrix(readdata('traffic_flux_transpose_row-normed', row.names=1))
    colnames(traffic_flux) <- areas
    rownames(traffic_flux) <- areas

    alt_traffic_flux <- array(0, dim=c(N,N,2))
    alt_traffic_flux[,,1] <- data.matrix(readdata('uk_forward_commute_flow', row.names=1))
    alt_traffic_flux[,,2] <- data.matrix(readdata('uk_reverse_commute_flow', row.names=1))
    colnames(alt_traffic_flux) <- areas
    rownames(alt_traffic_flux) <- areas

  })
  if (!is.null(env$opt$limit_area) && !is.null(env$opt$limit_radius)) {
    source("sandbox/limit_data.r")
    limit_data_by_distance(env, env$opt$limit_area, env$opt$limit_radius)
    print("Limited data to areas:")
    print(env$areas)
  }
  env
}