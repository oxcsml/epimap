##########################################################################
##########################################################################
# If run from shell script, process command line arguments and run
if (!interactive()) {
  source("dataprocessing/read_data.r")
  source("mapping/epimap.r")

  Rmap_opt = Rmap_options()
  Rmap_opt$cleaned_sample_id=1
  option_list = list(
    make_option(
      c("-s", "--spatialkernel"), 
      type="character", 
      default=Rmap_opt$spatialkernel,
      help=paste(
          "Use spatial kernel (matern12/matern32/matern52/exp_quad/none);",
          "default =", Rmap_opt$spatialkernel
      )
    ),
    make_option(
      c("-p", "--temporalkernel"), 
      type="character", 
      default=Rmap_opt$temporalkernel,
      help=paste(
          "Use temporal kernel (matern12/matern32/matern52/exp_quad/none);",
          "default =", Rmap_opt$temporalkernel
      )
    ),
    make_option(
      c("-l", "--localkernel"),
      type="character", 
      default=Rmap_opt$localkernel,
      help=paste("Use local kernel (local/none); default =", Rmap_opt$localkernel)
    ),
    make_option(
      c("-g", "--globalkernel"),  
      type="character",
      default=Rmap_opt$globalkernel,
      help=paste("Use global kernel (global/none); default =", Rmap_opt$globalkernel)
    ),
    make_option(
      c("-m", "--metapop"),
      type="character",
      default=Rmap_opt$metapop,
      help=paste(
          "metapopulation model for inter-region cross infections",
          "(none, or comma separated list containing radiation{1,2,3},traffic{forward,reverse},uniform,in,in_out);",
          "default = ", Rmap_opt$metapop
      )
    ),
    make_option(
      c("-o", "--observation"),
      type="character",
      default=Rmap_opt$observation,
      help=paste(
          "observation model",
          "(neg_binomial_{2,3}/poisson/cleaned_latent_sample/cleaned_latent_mean/cleaned_recon_sample);",
          "default =", Rmap_opt$observation
      )
    ),
    make_option(
      c("-x", "--cleaned_sample_id"),
      type="integer",
      default=Rmap_opt$cleaned_sample_id,
      help=paste("id of cleaned sample to use; default =", Rmap_opt$cleaned_sample_id)
    ),
    make_option(
      c("-c", "--chains"),
      type="integer",
      default=Rmap_opt$chains,
      help=paste("number of MCMC chains; default =", Rmap_opt$chains)
    ),
    make_option(
      c("-i", "--iterations"),
      type="integer",
      default=Rmap_opt$iterations,
      help=paste("Length of MCMC chains; defualt =", Rmap_opt$iterations)
    ),
    make_option(
      c("-n", "--time_steps"),
      type="integer",
      default=Rmap_opt$time_steps,
      help=paste("Number of periods to fit Rt in; default =",Rmap_opt$time_steps)
    ),
    make_option(
      c("-d", "--results_directory"),
      type="character",
      default=Rmap_opt$results_directory,
      help="If specified, store outputs in directory, otherwise use a unique directory"
    ),
    make_option(
      c("-r", "--clean_directory"),
      type="character",
      default=Rmap_opt$clean_directory,
      help="If specified, store outputs in directory, otherwise use a unique directory"
    ),
    make_option(
      c("-t", "--task_id"),
      type="integer",
      default=0,
      help="Task ID for Slurm usage. By default, turned off [0]."
    )
  ) 

  opt_parser = OptionParser(option_list=option_list)
  parsed_opt = parse_args(opt_parser)

  # If using Slurm, override other CLI options and use grid instead.
  if (parsed_opt$task_id > 0) {
    grid = expand.grid(
      spatialkernel=c("matern12", "matern32", "matern52", "exp_quad", "none"), 
      metapop=c("radiation1,uniform,in", "radiation1,uniform,in_out", "radiation2,uniform,in", "radiation2,uniform,in_out", "radiation3,uniform,in", "radiation3,uniform,in_out", "uniform,in", "uniform,in_out", "none"), 
      observation=c("neg_binomial_2", "neg_binomial_3", "poisson", "cleaned_sample","cleaned_mean"),
      localkernel=c("local","none"),
      globalkernel=c("global","none")
    )
    grid = sapply(grid, as.character)
    update = as.list(grid[parsed_opt$task_id, ]) 
    for (name in names(update)){
      parsed_opt[name] = update[name]
    }
  }

  for (o in names(parsed_opt)) {
    Rmap_opt[o] = parsed_opt[o]
  }

  env = Rmap_setup(Rmap_opt)
  
  Rmap_merge(env,c(1,2,3,4,5,6,7,8,9,10))
}
##########################################################################



