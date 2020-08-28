# Rmap

## Setup conda environment:
* `conda env create -f environment.yml`
* `conda activate Rmap`

## Current process:

In python3:
* `process_uk_cases.py` to produce `data/uk_cases.csv`
* `process_mobility.py` to produce `data/mobility.csv`
* `process_metadata.py` to produce `data/metadata.csv`

In R:
* `process_data.r` to produce `data/areas.csv` and `data/distances.csv`
* `run.r` to estimate Rt's, save to `fits/Rmap_[options]_[Rt/Cproj/logpred].csv`, `fits/Rmap_[options]_pairs.pdf` and ``fits/Rmap_[options]_stanfit.rds`
* File name convention: `Rmap_[datetime_id]_[options]` where `datetime_id` is a numeric string consisting of date and time when `run.r` is run, which hopefully serves as unique id.
* Options described in top of `run.r`
* (or, from Python, `run.py` to estimate Rt's and Ct projections)

From command line:
* `reinflate.sh [datetime_id]` To post-process data for website.

Visualisation
* Run `pip install --user -r requirements.txt` to install requirements
* Run `python3 run.py` to run the server
* Browse to `http://localhost:8000/?map=[datetime_id]`


Next steps to be tracked in Issues: https://github.com/rs-delve/Rmap/issues
