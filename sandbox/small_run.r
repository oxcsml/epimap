source('mapping/epimap.r')
source('sandbox/limit_data.r')

opt = Rmap_options(
    first_day_modelled = "2020-08-15",
    last_day_modelled = "2020-11-27",
    spatialkernel = "none",
    globalkernel = "none",
    days_ignored=NULL,
    days_per_step = 7,
    num_steps_forecasted = 21,
    iterations=3000,
    fixed_gp_time_length_scale = -1,
    fixed_gp_space_length_scale = -1,
    clean_directory='fits/clean-1205',
    results_directory='fits/test-daily-nb3-local'
)

env = Rmap_setup(opt)
limit_data_multi(env,list(
  list(area='Oxford',distance=.6),
  list(area='Birmingham',distance=.6)
))
#limit_data_by_distance(env,'Oxford',0.6)
# timing: 
# Oxford .8 -> 64 areas, ~ 10 hours
# Oxford .6 -> 26 areas, 2.8 hours
# Oxford .4 -> 8 areas, ~ .5 hours (removed local effects)
# Oxford .3 Birmingham .3 -> 21 areas, ~2.3 hours
# Oxford .5 Birmingham .5 -> 43 areas, ~ hours
# Oxford .6 Birmingham .6 -> 58 areas, ~ 9.2 hours

print(env$N)
print(env$areas)
print(env$N_region)
Rmap_run(env)
Rmap_postprocess(env)



# Rmap_merge(env,c(1,2,3))
