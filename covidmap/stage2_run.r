source('covidmap/stage2.r')
opt = covidmap_stage2_get_cmdline_options()
env = covidmap_stage2_setup(opt)
covidmap_stage2_run(env)
covidmap_stage2_postprocess(env)
