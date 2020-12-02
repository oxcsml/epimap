source('covidmap/stage1.r')

opt = covidmap_stage1_get_cmdline_options()
covidmap_stage1_combine(opt)

