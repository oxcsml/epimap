source('covidmap/stage2.r')
opt = covidmap_stage2_get_cmdline_options()
env = covidmap_stage2_setup(opt)
if (identical(env$opt$approximation,"twostage")) {
  singlearea_sample_ids = c(1,2,3,4,5,6,7,8,9,10)
  region_ids = c(0)
} else if (identical(env$opt$approximation,"regional")) {
  singlearea_sample_ids=c(0)
  region_ids=c(1,2,3,4,5,6,7,8,9)
}
covidmap_stage2_merge(env,
  singlearea_sample_ids=singlearea_sample_ids,
  region_ids=region_ids
)
