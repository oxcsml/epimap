source('mapping/epimap.r')
source('sandbox/limit_data.r')

opt = Rmap_options(
    iterations=3000,
    cleaned_sample_id = 1,
    observation_data='infections_latent',
    observation_model='gaussian'
)
env = Rmap_setup(opt)
limit_data_by_distance(env,'Manchester',.3)
Rmap_run(env)
Rmap_postprocess(env)


# Rmap_merge(env,c(1,2,3))
