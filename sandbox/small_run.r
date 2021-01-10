
system("cat sandbox/small_run.r")

source('mapping/epimap.r')
source('sandbox/limit_data.r')


opt = Rmap_options(
    stage = "full",
    region_id = 1,
    first_day_modelled = "2020-08-26",
    last_day_modelled = "2020-12-08",
    spatialkernel = "none",
    temporalkernel = "matern12",
    globalkernel = "none",
    localkernel = "local",
    days_ignored=NULL,
    days_per_step = 7,
    num_steps_forecasted = 3,
<<<<<<< HEAD
    iterations=30,
    fixed_gp_time_length_scale = -1,
    fixed_gp_space_length_scale = -1,
    clean_directory='fits/clean-1221',
    results_directory='fits/small-run'
=======
    iterations=3000,
    fixed_gp_time_length_scale = -1,
    fixed_gp_space_length_scale = -1,
    clean_directory='fits/clean-1215',
    results_directory='fits/test-weekly-northeast-newrt-local'
>>>>>>> 52e3216bcbca82f3c36a56d70a4d5b44214b47d9
)

env = Rmap_setup(opt)
#limit_data_multi(env,list(
  #list(area='Oxford',distance=.6),
  #list(area='Birmingham',distance=.6)
#  list(area='Westminster',distance=.8)
#))
#limit_data_by_distance(env,'Oxford',0.6)
# timing: 
# Oxford .8 -> 64 areas, ~ 10 hours
# Oxford .6 -> 26 areas, 2.8 hours
# Oxford .4 -> 8 areas, ~ .5 hours (removed local effects)
# Oxford .3 Birmingham .3 -> 21 areas, ~2.3 hours
# Oxford .5 Birmingham .5 -> 43 areas, ~ hours
# Oxford .6 Birmingham .6 -> 58 areas, ~ 9.2 hours

print(env$N)
print(env$areas[env$modelled_region[,1]==1])


Rmap_run(env)
Rmap_postprocess(env)



# Rmap_merge(env,c(1,2,3))
