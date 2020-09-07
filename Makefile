SHELL=/bin/bash
CONDAROOT=../miniconda3

environment:
	source $(CONDAROOT)/bin/activate && conda env create -f environment.yml
	source $(CONDAROOT)/bin/activate && conda activate Rmap && pip install -r requirements.txt
	git clone git@github.com:rs-delve/covid19_datasets.git
	source $(CONDAROOT)/bin/activate && conda activate Rmap && pip install -e covid19_datasets
		
preprocess-data:
	source $(CONDAROOT)/bin/activate && conda activate Rmap && python process_uk_cases.py
	source $(CONDAROOT)/bin/activate && conda activate Rmap && python process_mobility.py
	source $(CONDAROOT)/bin/activate && conda activate Rmap && python process_metadata.py
	source $(CONDAROOT)/bin/activate && conda activate Rmap && Rscript process_data.r
	# source $(CONDAROOT)/bin/activate && conda activate Rmap && Rscript process_radiation_fluxes.r
	source $(CONDAROOT)/bin/activate && conda activate Rmap && Rscript process_delays.r