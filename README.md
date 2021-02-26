# Rmap

## Installation

* Clone this repository with `git clone https://github.com/oxcsml/Rmap.git`.
* Install [conda](https://docs.conda.io/en/latest/) if not already installed.
* In the `Makefile`, set the `CONDAROOT` variable to the install path of conda. 
* Run `make environment`.

## Requirements

The requirements can be found in `environment.yml`, but should be automatically installed by the makefile/conda.

## Running the models:

First update the `CONDAROOT` variable in `Makefile` to point to the conda conda install.

Before running models, run `make preprocess-data` to pull the latest case data and regenerate various files. Note if the the file `data/uk_traffic.csv` already exists, then it will not be regenerated as this is very slow and does not need to be done regularly. To force a regeneration for whatever reason, simply delete the file.

NOTE: In the following, many of the scripts are designed for use with a SLURM job management system. If your system does not have one of these, then please consult the scripts for the required workflow, and adapt to your specific system. In general some kind of parallelism is required when running these models. Running each subtask sequentially would take prohibitively long.

To run a singlearea approximation, the script `scripts/submit-run-singlearea.sh` will automatically dispatch SLURM jobs to compute the models over the individual regions and then recombine them, saving the results and the cases file used in the specified directory, and using the options string specified.

After running a singlearea approximation, you can then run either a `twostage` or a `regional` approximation via the `scripts/submit-run-twostage.sh` and `scripts/submit-run-regional.sh` scripts, supplying the same results directory as for the single area approximation in order to load the correct singlearea results.

After running these models, the output director structure will look like:
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
and places the post processed results in `docs/assets/data`. To start the website locally see `docs/README.md`.

To run the models as per the website, see the script `scripts/daily_update.sh` for the workflow.

Next steps to be tracked in Issues: https://github.com/rs-delve/Rmap/issues

## Repo structure

* `covidmap` - contains code to run models specifically for COVID19 in the UK. The code here calls general models defined in `epimap` after loading the data and doing any preprocessing required for the specific scenario.
* `data` - contains data specific to the UK.
* `dataprocessing` - contains scripts for preprocessing specifc data about the UK.
* `doc` - contains misc documents about the project. No code.
* `docs` - contains the source code and data for the public website.
* `epimap` - contains the basic library of models used. These models could be applied more generally than to the UK, if supplied with the correct data.
* `evaluation` - contains code to evaluate the performance of models.
* `sandbox` - contains misc non-critical code.
* `scripts` - contains scripts for running models on UK data using SLURM job management.
