source('covidmap/stage1.r')

# Wrapper script to recombine the results produced
# inparallel from multiple bootstrap runs of stage1_run.r
opt = covidmap_stage1_get_cmdline_options()
covidmap_stage1_bootstrap_combine(opt)

