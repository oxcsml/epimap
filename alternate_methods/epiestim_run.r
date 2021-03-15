source("alternate_methods/epiestim.r")

# Wrapper script to run a single area under the 
# singlearea approximation
opt = epiestim_get_cmdline_options()
print(opt$area_index)
epiestim_run(opt$area_index, opt)

