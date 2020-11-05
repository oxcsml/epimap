
source('covidmap/stage1.r')
opt = covidmap_stage1_get_cmdline_options()
print(opt$area_index)
covidmap_stage1_run(opt$area_index, opt)

