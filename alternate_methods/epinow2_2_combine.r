source("alternate_methods/epinow2_2.r")

# Wrapper script to recombine the results produced
# inparallel from a run of stage1_run.r
opt = epinow2_get_cmdline_options()
epinow2_combine(opt)
