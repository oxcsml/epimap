source("dataprocessing/read_data.r")
source('mapping/epimap.r')
opt = epimap_get_cmdline_options()
env = Rmap_setup(opt)
if (identical(env$opt$stage,"map")) {
  cleaned_sample_ids = c(1,2,3,4,5,6,7,8,9,10)
  region_ids = c(0)
} else if (identical(env$opt$stage,"full")) {
  cleaned_sample_ids=c(0)
  region_ids=c(1,2,3,4,5,6,7,8,9)
}
Rmap_merge(env,
  cleaned_sample_ids=cleaned_sample_ids,
  region_ids=region_ids
)



