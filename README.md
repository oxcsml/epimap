# Rmap

Setup conda environment:
* conda env create -f environment.yml
* conda activate Rmap

Current process:
* process_uk_cases.py to produce uk_cases.csv
* process_mobility.py to produce mobility.csv
* process_metadata.py to produce metadata.csv
* run_cori.r to estimate Rt's, saved to RtCproj.csv

Visualisation
* Run `pip install --user -r requirements.txt` to install requirements
* Run `python app.py` to run the server
* Browse to `http://localhost:8050`

Next steps to be tracked in Issues: https://github.com/rs-delve/Rmap/issues
