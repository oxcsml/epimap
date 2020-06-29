# Rmap

Setup conda environment:
* conda env create -f environment.yml
* conda activate Rmap

Current process:
* process_uk_cases.py to produce uk_cases.csv
* run_cori.r to estimate Rt's, saved to Rt.csv

Visualisation
* Run `pip install --user -r requirements.txt` to install requirements
* Run `python app.py` to run the server
* Browse to `http://localhost:8050`

TODO Next steps:
* Code to plot nicely R maps.
* Obtain local authority coordinates, population density, population size, mobility data
* Get GP working
* Extend Cori et al method to account for delay in reporting, weekend effects
* Figure out how to incorporate both death and cases
* Figure out how to deal with more timely but messier data (e.g. NHS triage, Zoe app)
