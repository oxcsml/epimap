library("rjson")

#' Read in the data required to run the covid map models.
#' Loads basic data from opt$data_directory.
#' Loads cases data from opt$results_directory, ensuring that results and the case
#' counts used to create them are coupled.
#' Loads single area approximation results from opt$results_directory/singlearea in
#' order to run more complex models.
covidmap_read_data = function(env) { 
  with(env,{

    readdata = function(filename,...) {
      read.csv(sprintf('%s/%s.csv',opt$data_directory,filename),...)
    }
    readresults = function(filename,...) {
      read.csv(sprintf('%s/%s.csv',opt$results_directory,filename),...)
    }
    readsinglearea = function(filename,...) {
      read.csv(sprintf('%s/singlearea/%s.csv',opt$results_directory,filename),...)
    }

    #########################################################
    # Read infection profile data
    infprofile <- readdata("serial_interval")$fit
    Tip <- 30
    infprofile <- infprofile[1:Tip]
    infprofile <- infprofile/sum(infprofile)
    D <- length(infprofile)

    # Compute the testing delay profile
    Tdp <- 21
    presymptomdays <- 2
    Tdpnz <- Tdp - presymptomdays
    Adp <- 5.8
    Bdp <- 0.948
    testdelayprofile <- pgamma(1:Tdpnz,shape=Adp,rate=Bdp)
    testdelayprofile <- testdelayprofile/testdelayprofile[Tdpnz]
    testdelayprofile <- testdelayprofile - c(0.0,testdelayprofile[1:(Tdpnz-1)])
    testdelayprofile <- c(rep(0, presymptomdays), testdelayprofile)

    # Read non-temporal data about areas
    df <- readdata("areas",row.names=1)
    N <- nrow(df)
    areas <- rownames(df)
    quoted_areas <- sapply(areas,function(s) paste('"',s,'"',sep=''))
    geoloc <- df[,1:2]
    population <- df[,3]

    # NHS region data
    sparse_region <- df[,-3:-1, drop=FALSE]
    N_region <- ncol(sparse_region)
    region_df <- readdata("nhs_regions", row.names=1)
    regions <- rownames(region_df)
    quoted_regions <- sapply(regions,function(s) paste('"',s,'"',sep=''))

    # Read on precomputed distances between region centres.
    geodist <- readdata("distances",row.names=1)
    colnames(geodist) <- areas

    # Read case data from the results folder. This ensures that the case data
    # and results are generasted from are saved together
    cases <- readresults("cases")
    ind <- sapply(cases[,2], function(s)
        !(s %in% c('Outside Wales','Unknown','...17','...18'))
    )
    cases <- cases[ind,]
    AllCount <- cases[,3:ncol(cases)]
    Tall <- ncol(AllCount)
    dates <- as.Date(colnames(AllCount), format='X%Y.%m.%d')
    colnames(AllCount) <- dates
    rownames(AllCount) <- areas

    # All methods except the singlearea approximation require the results of a 
    # singlearea approximation to compute. If we are running such a method, load
    # a particular sample from the related singlearea results, run on the same cases
    # data.
    if (identical(opt$approximation, "regional") | identical(opt$approximation, "twostage")) {
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

      # Read data about how to split UK into regions to model for full model
      js = fromJSON(file=sprintf('%s/%s',opt$data_directory,"region-groupings.json")) ### TODO FIX
      inferred_region = sparse_region
      modelled_region = 0*sparse_region
      for (i in 1:opt$num_regions) {
        modelled_region[js[[i]],i]=1
      }
      stopifnot(all(inferred_region<=modelled_region))
    }

    # Load traffic flux data between regions
    traffic_flux <- array(0, dim=c(N,N,2))
    traffic_flux[,,1] <- data.matrix(readdata('traffic_flux_row-normed', row.names=1))
    traffic_flux[,,2] <- data.matrix(readdata('traffic_flux_transpose_row-normed', row.names=1))
    colnames(traffic_flux) <- areas
    rownames(traffic_flux) <- areas

    # Load alternative commuter flow data for intra region flux
    alt_traffic_flux <- array(0, dim=c(N,N,2))
    alt_traffic_flux[,,1] <- data.matrix(readdata('uk_forward_commute_flow', row.names=1))
    alt_traffic_flux[,,2] <- data.matrix(readdata('uk_reverse_commute_flow', row.names=1))
    colnames(alt_traffic_flux) <- areas
    rownames(alt_traffic_flux) <- areas

  })
  # If desired, limit the number of regions used in the computation based on a central region
  # and a radius around it.
  if (!is.null(env$opt$limit_area) && !is.null(env$opt$limit_radius)) {
    source("sandbox/limit_data.r")
    limit_data_by_distance(env, env$opt$limit_area, env$opt$limit_radius)
    print("Limited data to areas:")
    print(env$areas)
  }
  env
}
