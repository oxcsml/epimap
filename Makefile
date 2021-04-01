SHELL=/bin/bash
CONDAROOT=/data/ziz/not-backed-up/bhe/miniconda3

environment:
	source $(CONDAROOT)/bin/activate && conda env create -f environment.yml
	git clone git@github.com:rs-delve/covid19_datasets.git dataprocessing
	source $(CONDAROOT)/bin/activate && conda activate Rmap && pip install -e dataprocessing/covid19_datasets
		
preprocess-data:
	cd dataprocessing/covid19_datasets && git pull && cd -
	source $(CONDAROOT)/bin/activate && conda activate Rmap && python dataprocessing/process_uk_cases.py
	source $(CONDAROOT)/bin/activate && conda activate Rmap && python dataprocessing/process_mobility.py
	source $(CONDAROOT)/bin/activate && conda activate Rmap && python dataprocessing/process_metadata.py
	source $(CONDAROOT)/bin/activate && conda activate Rmap && python dataprocessing/process_site_data.py
	source $(CONDAROOT)/bin/activate && conda activate Rmap && python dataprocessing/process_region_site_data.py
	source $(CONDAROOT)/bin/activate && conda activate Rmap && python scripts/group_local_authorities.py data/traffic_flux_row-normed.csv data/traffic_flux_transpose_row-normed.csv data/areas.csv data/region-groupings.json --threshold=0.8
	source $(CONDAROOT)/bin/activate && conda activate Rmap && Rscript dataprocessing/process_data.r
	# source $(CONDAROOT)/bin/activate && conda activate Rmap && Rscript process_radiation_fluxes.r
	source $(CONDAROOT)/bin/activate && conda activate Rmap && python dataprocessing/traffic/create_traffic_matrix.py
	source $(CONDAROOT)/bin/activate && conda activate Rmap && python dataprocessing/traffic/process_alt_traffic.py
