library(EpiNow2)
library(data.table)
library(optparse)


option_list <- list(
            make_option(
                    "--first_day_modelled",
                    action = "store",
                    default = NA,
                    type = "character",
                    help = "First date of date to use"
                ),
            make_option(
                    "--weeks_modelled",
                    action = "store",
                    default = NA,
                    type = "integer",
                    help = "Number of weeks to model"
                ),
            make_option(
                    "--forecast_horizon",
                    action = "store",
                    default = NA,
                    type = "integer",
                    help = "Number of days for which to forecast"
                ),
            make_option(
                    "--area",
                    action = "store",
                    default = NA,
                    type = "character",
                    help = "Name of the area to model, case sensitive!"
                ),
            make_option(
                    "--case_counts",
                    action = "store",
                    default = NA,
                    type = "character",
                    help = "Path to preprocessed input data"
                ),
            make_option(
                    "--ncores",
                    action = "store",
                    default = NA,
                    type = "integer",
                    help = "Number of cores to tell stan to use"
                ),
            make_option(
                    "--output_folder",
                    action = "store",
                    default = NA,
                    type = "character",
                    help = "Where to store raw results"
                ),
            make_option(
                    "--region_codes",
                    action = "store",
                    default = NA,
                    type = "character",
                    help = "Path to JSON file with region codes. Format area_name -> code."
                ),
            make_option(
                    "--verbose",
                    action = "store_true",
                    default = FALSE,
                    type = "logical",
                    help = "Whether stan should be verbose"
                ),
            make_option(
                    "--area_index",
                    action="store",
                    type="integer",
                    default=0,
                    help="Area index; used if area not provided."
                )
        )

opt <- parse_args(OptionParser(option_list = option_list))

start_date <- as.IDate(opt$first_day_modelled) 
df <- fread(opt$case_counts)
region_codes <- rjson::fromJSON(file = opt$region_codes)
if(is.na(opt$area)){
    if (opt$area_index==0) {
        stop("Area index 0.")
    }
    opt$area = names(region_codes)[opt$area_index]
}
area_code <- region_codes[[opt$area]]
end_date <- start_date + opt$weeks_modelled * 7
region_data <- subset(df[df$region == opt$area], select = c("date", "confirm"))
input_data <- region_data[region_data$date >= start_date & region_data$date < end_date] 
input_data$date <- as.Date(input_data$date) # for compat with epinow2
print(tail(input_data, 1))

# as close as we can reasonably get to the paper
reporting_delay <- list(log(6.5^2 / sqrt(17^2 + 6.5^2)),
                                0,
                                log(17^2 / 6.5^2 + 1),
                                    0,
                                    30)
names(reporting_delay) <- c("mean", "mean_sd", "sd", "sd_sd", "max")

# generation interval, incubation time and rt prior are as in the paper
generation_time <- get_generation_time(disease = "SARS-CoV-2", source = "ganyani")
incubation_period <- get_incubation_period(disease = "SARS-CoV-2", source = "lauer")
rt_options <- rt_opts(prior=list(mean=1, sd=1))

stan_options <- stan_opts(cores = opt$ncores)

# the setting of zero_threshold is to stop backfilling
estimates <- epinow(reported_cases = input_data,
                    generation_time = generation_time,
                    delays = delay_opts(incubation_period, reporting_delay),
                    rt = rt_options,
                    stan = stan_options,
                    horizon = opt$forecast_horizon,
                    verbose = opt$verbose,
                    zero_threshold = Inf 
                    )

standata <- rstan::extract(estimates$estimates$fit)
fitted_r_samples <- standata$R 
fitted_case_samples <- standata$imputed_reports

write.table(
            fitted_r_samples,
            file = file.path(opt$output_folder, sprintf("%d_r_samples.txt", area_code)),
            row.names = FALSE,
            col.names = FALSE
        )

write.table(
            fitted_case_samples,
            file = file.path(opt$output_folder, sprintf("%d_case_samples.txt", area_code)),
            row.names = FALSE,
            col.names = FALSE
        )
