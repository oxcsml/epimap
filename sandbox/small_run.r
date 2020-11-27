source('mapping/epimap.r')
source('sandbox/limit_data.r')

opt = Rmap_options(
    iterations=3000,
    cleaned_sample_id = 2,
    observation_data='cleaned_latent_sample',
    observation_model='gaussian',
    fixed_gp_time_length_scale = 100.0,
    fixed_gp_space_length_scale = 3.0,
    clean_directory='fits/clean-1109',
    results_directory='fits/test-stochastic-forward-sampling'
)
env = Rmap_setup(opt)
limit_data_by_distance(env,'Manchester',.5)
Rmap_run(env)
Rmap_postprocess(env)


# Rmap_merge(env,c(1,2,3))
