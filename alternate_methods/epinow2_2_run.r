source("alternate_methods/epinow2_2.r")

# Wrapper script to run a single area under the 
# epinow2 model
opt = epinow2_get_cmdline_options()
print(opt$area_index)
epinow2_run(opt$area_index, opt)

