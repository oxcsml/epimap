# Rmap

Repository containing the code for the models used by the website [localcovid.info](https://localcovid.info).

## Repo structure

* `epimap` - contains the basic library of models used. These models could be applied more generally than to the UK, if supplied with the correct data.
* `covidmap` - contains code to run models specifically for COVID19 in the UK. The code here calls general models defined in `epimap` after loading the data and doing any preprocessing required for the specific scenario.
* `data` - contains data specific to the UK.
* `dataprocessing` - contains scripts for preprocessing specifc data about the UK.
* `doc` - contains misc documents about the project. No code.
* `docs` - contains the source code and data for the public website.
* `evaluation` - contains code to evaluate the performance of models.
* `sandbox` - contains misc non-critical code.
* `scripts` - contains scripts for running models on UK data using SLURM job management.


Next steps to be tracked in Issues: https://github.com/rs-delve/Rmap/issues

## Models

A brief description of the models are listed here. For more info please see the technical report. 

### Singlearea

Runs a renewal equation model on an individual area, using an AR1 process / Matern12 kernel Gaussian Process as the prior for the Rt through time.

Dispatch via `slurm/submit-run-singlearea.sh`

For parameter options, see `covidmap/stage1.r:covidmap_stage1_options`

### Twostage

Using a samples from the inferred underlying infections in each area (computed by the singlearea approximaton), this model inferes the Rt from these underlying infections in each area. A Kronecker factored Gaussian Process prior is placed over the Rt. It also incorporated metapopulation effects to include intra area transfer of cases. This effect is based on commuter flow data, and the prevalence of this intra area transfer is inferred. 

We infer over 10 samples of the underlying infections and then recombine the samples from each inference.

Dispatch via `slurm/submit-run-twostage.sh`, pointing to the same results directory as the single area results.

For parameter options, see `covidmap/stage2.r:covidmap_stage2_options`


### Regional

This model breaks the country into a number of regions for more efficient computation. It computes the Rt for areas in the specified region from the cases observed in the areas. In addition to performing inference on the areas in the region, additional areas are modelled to get a more accurate approximation. The additional areas added are the areas the contribute the top 80% of commuter flow into the region. The infered values of Rt for these additionally modelled regions are thrown away at the end. The metapopulation effects for areas not included in this extra modelling set are approximated from the results of a singlearea approximation of each of those areas.

Dispatch via `slurm/submit-run-regional.sh`, pointing to the same results directory as the single area results.

For parameter options, see `covidmap/stage2.r:covidmap_stage2_options`

## Installation

* Clone this repository with `git clone https://github.com/oxcsml/Rmap.git`.
* Install [conda](https://docs.conda.io/en/latest/) if not already installed.
* In the `Makefile`, set the `CONDAROOT` variable to the install path of conda. 
* Run `make environment`.


If you wish to run the alternative methods, EpiEstim must be installed manually by running `install.packages("EpiEstim")` in R.

## Requirements

The requirements can be found in `environment.yml`, but should be automatically installed by the makefile/conda.

## Running the models:

First update the `CONDAROOT` variable in `Makefile` to point to the conda conda install.

Before running models, run `make preprocess-data` to pull the latest case data and regenerate various files. Note if the the file `data/uk_traffic.csv` already exists, then it will not be regenerated as this is very slow and does not need to be done regularly. To force a regeneration for whatever reason, simply delete the file.

NOTE: In the following, many of the scripts are designed for use with a SLURM job management system. If your system does not have one of these, then please consult the scripts for the required workflow, and adapt to your specific system. In general some kind of parallelism is required when running these models. Running each subtask sequentially would take prohibitively long.

To run a singlearea approximation, the script `scripts/submit-run-singlearea.sh` will automatically dispatch SLURM jobs to compute the models over the individual regions and then recombine them, saving the results and the cases file used in the specified directory, and using the options string specified.

After running a singlearea approximation, you can then run either a `twostage` or a `regional` approximation via the `scripts/submit-run-twostage.sh` and `scripts/submit-run-regional.sh` scripts, supplying the same results directory as for the single area approximation in order to load the correct singlearea results.

After running these models, the output directory structure will look like:
```
+-- singlearea
|   +-- ...
+-- regional
|   +-- ...
+-- twostage
|   +-- ...
+-- cases.csv
```

After running the models, the data can be postprocessed to be displayed on the webview with the script `dataprocessing/reinflate.sh`. The script takes arguments
```
dataprocessing/reinflate.sh [path to results output] [name of postprocessed results]
```
and places the post processed results in `docs/assets/data`. In order to update the underlying cases data on the website run `python dataprocessing/process_site_data.py`. To start the website locally see `docs/README.md`.

To run the models as per the website, see the script `scripts/daily_update.sh` for the workflow.

