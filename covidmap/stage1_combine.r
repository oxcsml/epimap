source('covidmap/stage1.r')

# Wrapper script to recombine the results produced
# inparallel from a run of stage1_run.r
opt = covidmap_stage1_get_cmdline_options()
covidmap_stage1_combine(opt)

