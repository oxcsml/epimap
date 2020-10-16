source("dataprocessing/read_data.r")
source('mapping/epimap.r')
opt = epimap_get_cmdline_options()
env = Rmap_setup(opt)
Rmap_run(env)
Rmap_postprocess(env)


