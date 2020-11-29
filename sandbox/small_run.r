source('mapping/epimap.r')
source('sandbox/limit_data.r')

opt = Rmap_options(
    iterations=4000,
    cleaned_sample_id = 1,
    observation_data='latent_reports', #'cleaned_latent_sample',
    observation_model='negative_binomial_2', #'gaussian',
    fixed_gp_time_length_scale = 100.0,
    fixed_gp_space_length_scale = 3.0,
    clean_directory='fits/clean-1109',
    results_directory='fits/test-smallrun'
)
env = Rmap_setup(opt)
limit_data_by_distance(env,'Oxford',.8)
Rmap_run(env)
Rmap_postprocess(env)


# Rmap_merge(env,c(1,2,3))
