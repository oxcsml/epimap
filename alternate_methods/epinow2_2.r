library(EpiNow2)
library(data.table)
library(optparse)
library(rstan)
source("covidmap/utils.r")

epinow2_options = function(
    first_day_modelled   = NULL,
    last_day_modelled    = NULL,
    weeks_modelled       = NULL,
    days_ignored         = NULL,
    days_per_step        = 7,
    days_predicted       = 0,
    num_steps_forecasted = 3,

    num_iterations      = 2000,
    num_chains          = 4,
    data_directory      = 'data',
    results_directory   = 'fits/test',
    verbose             = TRUE,

    area_index          = 0
) {
  as.list(environment())
}

epinow2_run = function(area_index, opt) {
    print(opt)

    cases = read.csv(sprintf('%s/cases.csv',opt$results_directory))
    areas = rownames(read.csv(sprintf('%s/areas.csv',opt$data_directory),row.names=1))
    area = areas[area_index]
    dates = as.Date(names(cases)[3:ncol(cases)], format="X%Y.%m.%d")

    list[Nstep, Tstep, Tcond, Tlik, Tcur, Tignore] = process_dates_modelled(
        dates, 
        first_day_modelled = opt$first_day_modelled,
        last_day_modelled  = opt$last_day_modelled,
        days_ignored       = opt$days_ignored,
        weeks_modelled     = opt$weeks_modelled,
        days_per_step      = opt$days_per_step
    )

    start_date <- dates[Tcond + 1] 
    end_date <- dates[Tcond + (Nstep * opt$days_per_step)]
    days_forecast = Tstep * opt$num_steps_forecasted

    input_data = data.frame(date = dates, confirm = as.numeric(cases[area_index, 3:ncol(cases)]))
    input_data = input_data[input_data$date >= start_date & input_data$date <= end_date,]

    # These are just the default parameters that the package suggests in the README
    reporting_delay <- estimate_delay(rlnorm(1000,  log(3), 1), max_value = 15, bootstraps = 1)
    generation_time <- get_generation_time(disease = "SARS-CoV-2", source = "ganyani")
    incubation_period <- get_incubation_period(disease = "SARS-CoV-2", source = "lauer")
    rt_options <- rt_opts(prior = list(mean = 2, sd = 0.2))
    stan_options <- stan_opts(cores = opt$num_chains, samples=opt$num_iterations)
    # setup_logging(threshold = "DEBUG", mirror_to_console = TRUE, file = file.path(opt$results_directory, "epinow2", "output", sprintf("log_%s.txt", area)), name = "EpiNow2")
    # setup_logging(threshold = "DEBUG", mirror_to_console = TRUE, file = file.path(opt$results_directory, "epinow2", "output", sprintf("log_%s.txt", area)), name = "Epinow2.epinow")
    # setup_logging(threshold = "DEBUG", mirror_to_console = TRUE, file = file.path(opt$results_directory, "epinow2", "output", sprintf("log_%s.txt", area)), name = "EpiNow2.epinow.estimate_infections")
    # setup_logging(threshold = "DEBUG", mirror_to_console = TRUE, file = file.path(opt$results_directory, "epinow2", "output", sprintf("log_%s.txt", area)), name = "EpiNow2.epinow.estimate_infections.fit")
    # setup_logging(threshold = "DEBUG", mirror_to_console = TRUE, file = file.path(opt$results_directory, "epinow2", "output", sprintf("log_%s.txt", area)), name = "EpiNow2")
    
    estimates <- epinow(reported_cases = input_data,
                        generation_time = generation_time,
                        delays = delay_opts(incubation_period, reporting_delay),
                        rt = rt_options,
                        stan = stan_options,
                        horizon = days_forecast,
                        verbose = opt$verbose)

    standata <- rstan::extract(estimates$estimates$fit)
    fitted_r_samples <- standata$R 
    fitted_case_samples <- standata$imputed_reports

    write.table(
            fitted_r_samples,
            file = file.path(opt$results_directory, "epinow2", "samples", sprintf("%s_r_samples.txt", area)),
            row.names = FALSE,
            col.names = FALSE
        )

    write.table(
                fitted_case_samples,
                file = file.path(opt$results_directory, "epinow2", "samples", sprintf("%s_case_samples.txt", area)),
                row.names = FALSE,
                col.names = FALSE
            )

    saveRDS(estimates$estimates$fit, file.path(opt$results_directory, "epinow2", "stanfits", sprintf("%s.rds", area)))
}

epinow2_combine = function(opt) {

    print(opt)

    area_date_dataframe <- function(areas,dates,provenance,data,data_names) {
        numareas <- length(areas)
        numdates <- length(dates)
        dates <- rep(dates,numareas)
        dim(dates) <- c(numareas*numdates)
        provenance <- rep(provenance,numareas)
        dim(provenance) <- c(numareas*numdates)
        areas <- rep(areas,numdates)
        dim(areas) <- c(numareas,numdates)
        areas <- t(areas)
        dim(areas) <- c(numareas*numdates)
        df <- data.frame(area=areas,Date=dates,data=data,provenance=provenance)
        colnames(df)[3:(ncol(df)-1)] <- data_names
        df
    }

    cases = read.csv(sprintf('%s/cases.csv',opt$results_directory))
    areas = rownames(read.csv(sprintf('%s/areas.csv',opt$data_directory),row.names=1))
    quoted_areas <- sapply(areas,function(s) paste('"',s,'"',sep=''))
    dates = as.Date(names(cases)[3:ncol(cases)], format="X%Y.%m.%d")

    list[Nstep, Tstep, Tcond, Tlik, Tcur, Tignore] = process_dates_modelled(
        dates, 
        first_day_modelled = opt$first_day_modelled,
        last_day_modelled  = opt$last_day_modelled,
        days_ignored       = opt$days_ignored,
        weeks_modelled     = opt$weeks_modelled,
        days_per_step      = opt$days_per_step
    )

    start_date <- dates[Tcond + 1] 
    days_modelled = Nstep * Tstep
    end_date <- dates[Tcond + days_modelled]
    days_forecast = Tstep * opt$num_steps_forecasted
    dates_modelled = dates[dates >= start_date & dates <= end_date]
    dates_predicted = seq(end_date,by=1,length.out=days_forecast)
    dates = c(dates_modelled, dates_predicted)

    provenance <- c(rep('inferred',days_modelled),rep('projected',days_forecast))

    N = length(areas)
    Nsample = 1
    
    ###
    # N = 10
    # areas = areas[1:N]
    # quoted_areas = quoted_areas[1:N]

    
    percentiles = c(.025,.1,.2,.25,.3,.4,.5,.6,.7,.75,.8,.9,.975)
    str_percentiles = c("2.5%","10%","20%","25%","30%","40%","50%","60%","70%","75%","80%","90%","97.5%")
    num_percentiles = length(percentiles)
    num_samples = opt$num_iterations

    Rt_percentiles = array(NA, c(N, days_modelled + days_forecast, num_percentiles))
    Cpred = array(NA, c(N, days_modelled, num_percentiles))
    Cproj = array(NA, c(N, days_forecast, num_percentiles))

    Xpred = array(NA, c(N, days_modelled, num_percentiles))
    Xproj = array(NA, c(N, days_forecast, num_percentiles))

    for (area_index in 1:N) {
        area <- areas[area_index]
        print(area)
    
        fit <- readRDS(paste(opt$results_directory, '/epinow2/stanfits/',area,'.rds',sep=''))
    
        skip <- num_samples / 2 / Nsample
        area_rt = summary(fit, pars = "R", probs=percentiles)$summary
        reports = summary(fit, pars = "reports", probs=percentiles)$summary
        infections = summary(fit, pars = "infections", probs=percentiles)$summary
        # Need to deal with the fact that epinow2 sometimes cuts off days at the start before the first case appears.
        total_days = days_modelled + days_forecast
        missing_days = total_days - dim(area_rt)[1]
        seeding_days = dim(infections)[1] - dim(reports)[1]
        # message("area_rt dims = ",dim(area_rt))
        for (p in 1:num_percentiles) {
            Rt_percentiles[area_index,(missing_days+1):total_days,p] = area_rt[,str_percentiles[p]]
        }

        Cpred[area_index,(missing_days+1):days_modelled,] = reports[1:(days_modelled - missing_days),str_percentiles]
        Cproj[area_index,,] = reports[(days_modelled - missing_days + 1):(total_days - missing_days),str_percentiles]
        Xpred[area_index,(missing_days+1):days_modelled,] = infections[(seeding_days + 1):(seeding_days + days_modelled - missing_days),str_percentiles]
        Xproj[area_index,1:days_forecast,] = infections[(seeding_days + days_modelled - missing_days + 1):(seeding_days + total_days - missing_days),str_percentiles]
           
    }

    Rt = aperm(Rt_percentiles, c(2,1,3))
    dim(Rt) <- c(N*(days_modelled + days_forecast), num_percentiles)
    Rt <- area_date_dataframe(
        quoted_areas,
        dates,
        provenance,
        format(round(Rt,2),nsmall=2),
        c("Rt_2_5","Rt_10","Rt_20","Rt_25","Rt_30","Rt_40","Rt_50",
        "Rt_60","Rt_70","Rt_75","Rt_80","Rt_90","Rt_97_5")
    )
    write.csv(Rt, paste(opt$results_directory, "/epinow2/Rt.csv", sep=""), quote=FALSE, row.names=FALSE)

    Cpred = aperm(Cpred, c(2,1,3))
    dim(Cpred) <- c(N*days_modelled, num_percentiles)
    Cpred <- area_date_dataframe(
        quoted_areas,
        dates_modelled,
        rep('inferred',days_modelled),
        Cpred,
        c("C_025","C_10","C_20","C_25","C_30","C_40","C_50",
            "C_60", "C_70","C_75","C_80","C_90","C_975")
    )
    write.csv(Cpred, paste(opt$results_directory, "/epinow2/Cpred.csv", sep=""), quote=FALSE, row.names=FALSE)


    Cproj = aperm(Cproj, c(2,1,3))
    dim(Cproj) <- c(N*days_forecast, num_percentiles)
    Cproj <- area_date_dataframe(
        quoted_areas,
        dates_predicted,
        rep('projected',days_forecast),
        Cproj,
        c("C_025","C_10","C_20","C_25","C_30","C_40","C_50",
            "C_60", "C_70","C_75","C_80","C_90","C_975")
    )
    write.csv(Cproj, paste(opt$results_directory, "/epinow2/Cproj.csv", sep=""), quote=FALSE, row.names=FALSE)

    Xpred = aperm(Xpred, c(2,1,3))
    dim(Xpred) <- c(N*days_modelled, num_percentiles)
    Xpred <- area_date_dataframe(
        quoted_areas,
        dates_modelled,
        rep('inferred', days_modelled),
        Xpred,
        c("X_025","X_10","X_20","X_25","X_30","X_40","X_50",
            "X_60", "X_70","X_75","X_80","X_90","X_975")
    )
    write.csv(Xpred, paste(opt$results_directory, "/epinow2/Xpred.csv", sep=""), quote=FALSE, row.names=FALSE)


    Xproj = aperm(Xproj, c(2,1,3))
    dim(Xproj) <- c(N*days_forecast, num_percentiles)
    Xproj <- area_date_dataframe(
        quoted_areas,
        dates_predicted,
        rep('projected',days_forecast),
        Xproj,
        c("X_025","X_10","X_20","X_25","X_30","X_40","X_50",
            "X_60", "X_70","X_75","X_80","X_90","X_975")
    )
    write.csv(Xproj, paste(opt$results_directory, "/epinow2/Xproj.csv", sep=""), quote=FALSE, row.names=FALSE)



}

epinow2_cmdline_options = function(opt = epinow2_options()) {
    list(
    make_option(
        c("--first_day_modelled"),
        type="character",
        default=opt$first_day_modelled,
        help=paste("Date of first day to model; default =",opt$first_day_modelled)
    ),
    make_option(
        c("--weeks_modelled"),
        type="integer",
        default=opt$weeks_modelled,
        help=paste("Number of weeks to model; default =",opt$weeks_modelled)
    ),
    make_option(
        c("--last_day_modelled"),
        type="character",
        default=opt$last_day_modelled,
        help=paste("Date of last day to model; default =",opt$last_day_modelled)
    ),
    make_option(
        c("--days_ignored"),
        type="integer",
        default=opt$days_ignored,
        help=paste("Number of recent days ignored, default",opt$days_ignored)
    ),
    make_option(
        c("--days_per_step"),
        type="integer",
        default=opt$days_per_step,
        help=paste("Days per modelling step, default",opt$days_per_step)
    ),
    make_option(
        c("--days_predicted"),
        type="integer",
        default=opt$days_predicted,
        help=paste("Days predicted; default =",opt$days_predicted)
    ),
    make_option(
        c("--num_steps_forecasted"),
        type="integer",
        default=opt$num_steps_forecasted,
        help=paste("Number of steps to forecast, default",opt$num_steps_forecasted)
    ),
    make_option(
        "--area_index",
        action = "store",
        default = opt$area_index,
        type = "integer",
        help = "Number of the area to model in the cases file"
    ),
    make_option(
        "--num_chains",
        action = "store",
        default = opt$num_chains,
        type = "integer",
        help = "Number of cores to tell stan to use"
    ),
    make_option(
        "--num_iterations",
        action = "store",
        default = opt$num_iterations,
        type = "integer",
        help = "Number of cores to tell stan to use"
    ),
    make_option(
        c("--results_directory"), 
        type="character", 
        default=opt$results_directory, 
        help=paste("Directory to put cleaned results in, default ",opt$results_directory)
    ),
    make_option(
        c("--data_directory"), 
        type="character", 
        default=opt$data_directory, 
        help=paste("Directory to get data from, default ",opt$data_directory)
    ),
    make_option(
        "--verbose",
        action = "store_true",
        default = FALSE,
        type = "logical",
        help = "Whether stan should be verbose"
    )
    )
}

epinow2_get_cmdline_options = function(opt=epinow2_options()) {
  cmdline_opt = epinow2_cmdline_options(opt)
  opt_parser = OptionParser(option_list=cmdline_opt)
  parsed_opt = parse_args(opt_parser)
  for (o in names(parsed_opt)) {
    opt[o] = parsed_opt[o]
  }
  opt
}