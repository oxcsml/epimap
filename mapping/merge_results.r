source("dataprocessing/read_data.r")
source('mapping/epimap.r')
opt = epimap_get_cmdline_options()
env = Rmap_setup(opt)
Rmap_merge(env,c(1,2,3,4,5,6,7,8,9,10))



