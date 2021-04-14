sample_directory=simulation/latent_epidemic/tehtropolis/sample_rss_paper_1
results_directory=simulation_fits/sample_rss_paper_1

mkdir -p $results_directory
cp $sample_directory/cases.csv $results_directory

common_options="\
    --data_directory $sample_directory \
    --first_day_modelled 2020-04-09 \
    --days_ignored 21 \
    --days_predicted 21 \
    --num_steps_forecasted 3"

epimap_stage1_options="\
    --fixed_gp_time_length_scale 200"

epimap_stage2_options="\
    --fixed_gp_space_length_scale 0.05 \
    --fixed_gp_time_length_scale 200"

./slurm/submit-run-epiestim.sh \
    $results_directory \
    "$common_options" \
    5

./slurm/submit-run-epinow2.sh \
    $results_directory \
    "$common_options" \
    5

./slurm/submit-run-singlearea.sh \
    $results_directory \
    "$common_options $epimap_stage1_options" \
    5

./slurm/submit-run-twostage.sh \
    $results_directory \
    "$common_options $epimap_stage2_options" \
    10 

./slurm/submit-run-regional.sh \
    $results_directory \
    "$common_options $epimap_stage2_options" \
    1



library(EpiEstim) 
library(optparse) 
library(gsubfn) 
library(plyr) 
library(data.table)  
source("covidmap/read_data.r") 
source("covidmap/utils.r")   
source("alternate_methods/epiestim.r")                                                                                                                                       

opt = epiestim_options()                                                                                                                                                     
opt$data_directory = "simulation/latent_epidemic/tehtropolis/sample_rss_paper_1"     
opt$results_directory = "simulation_fits/sample_rss_paper_1"                                                                                                    
opt$first_day_modelled = "2020-04-09"                                                                                                                                        
opt$days_ignored = 21                                                                                                                                                        
opt$days_predicted = 21                                                                                                                                                      
opt$num_steps_forecast = 3    
opt$num_samples = 10 

# opt = epinow2_options()
# opt$data_directory = "simulation/latent_epidemic/tehtropolis/sample"     
# opt$results_directory = "simulation_fits/test"                                                                                                    
# opt$first_day_modelled = "2020-04-09"                                                                                                                                        
# opt$days_ignored = 21                                                                                                                                                        
# opt$days_predicted = 21                                                                                                                                                      
# opt$num_steps_forecast = 3    
# opt$num_samples = 10 

# source('covidmap/stage2.r')
# opt = covidmap_stage2_get_cmdline_options()
# opt$data_directory = "simulation/latent_epidemic/tehtropolis/sample"     
# opt$results_directory = "simulation_fits/test"                                                                                                    
# opt$first_day_modelled = "2020-04-09"                                                                                                                                        
# opt$days_ignored = 21                                                                                                                                                        
# opt$days_predicted = 21                                                                                                                                                      
# opt$num_steps_forecast = 3    
# opt$num_samples = 10 
# opt$approximation="twostage"
# opt$singlearea_sample_id=1
# opt$region_id=0

