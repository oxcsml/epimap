source('run.r')
source('limit_data.r')

for (i in 1:3) {
  opt = Rmap_options(iterations=60,cleaned_sample_id=i)
  env = Rmap_setup(opt)
  limit_data_by_distance(env,'Manchester',.3)
  Rmap_run(env)
  Rmap_postprocess(env)
}

Rmap_merge(env,c(1,2,3))
