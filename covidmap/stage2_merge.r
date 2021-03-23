source('covidmap/stage2.r')
opt = covidmap_stage2_get_cmdline_options()
env = covidmap_stage2_setup(opt)
if (identical(env$opt$approximation,"twostage")) {
  singlearea_sample_ids = 1:env$opt$num_samples
  region_ids = c(0)
} else if (identical(env$opt$approximation,"regional")) {
  singlearea_sample_ids=c(0)
  region_ids=1:env$opt$num_regions
}
covidmap_stage2_merge(env,
  singlearea_sample_ids=singlearea_sample_ids,
  region_ids=region_ids
)
