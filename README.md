# Rmap

Setup conda environment:
* `conda env create -f environment.yml`
* `conda activate Rmap`

Current process:
* `process_uk_cases.py` to produce `data/uk_cases.csv`
* `process_mobility.py` to produce `data/mobility.csv`
* `process_metadata.py` to produce `data/metadata.csv`
* `process_data.r` to produce `data/areas.csv` and `data/distances.csv`
* `run.r` to estimate Rt's, save to `fits/RtCproj_[options].csv` and `fits/stanfit_[options].rds`
* (or, from Python, `run.py` to estimate Rt's and Ct projections)

Visualisation
* Run `pip install --user -r requirements.txt` to install requirements
* Run `python3 run.py` to run the server
* Browse to `http://localhost:8000`

Next steps to be tracked in Issues: https://github.com/rs-delve/Rmap/issues
