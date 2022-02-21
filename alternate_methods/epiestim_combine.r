source("alternate_methods/epiestim.r")

# Wrapper script to recombine the results produced
# inparallel from a run of stage1_run.r
opt = epiestim_get_cmdline_options()
epiestim_combine(opt)

