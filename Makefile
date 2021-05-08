SHELL=/bin/bash
CONDAROOT=/data/ziz/not-backed-up/teh/miniconda3

environment:
	source $(CONDAROOT)/bin/activate && conda env create -f environment.yml
	#git clone git@github.com:epimap/covid19_datasets.git dataprocessing/covid19_datasets
		
preprocess-data:
	cd dataprocessing/covid19_datasets && git pull && cd -
	source $(CONDAROOT)/bin/activate && conda activate Rmap-daily-update && python dataprocessing/process_uk_cases.py
	source $(CONDAROOT)/bin/activate && conda activate Rmap-daily-update && python dataprocessing/process_mobility.py
	source $(CONDAROOT)/bin/activate && conda activate Rmap-daily-update && python dataprocessing/process_metadata.py
	source $(CONDAROOT)/bin/activate && conda activate Rmap-daily-update && python dataprocessing/process_site_data.py
	source $(CONDAROOT)/bin/activate && conda activate Rmap-daily-update && python dataprocessing/process_region_site_data.py
	source $(CONDAROOT)/bin/activate && conda activate Rmap-daily-update && python scripts/group_local_authorities.py data/traffic_flux_row-normed.csv data/traffic_flux_transpose_row-normed.csv data/areas.csv data/region-groupings.json --threshold=0.8
	source $(CONDAROOT)/bin/activate && conda activate Rmap-daily-update && Rscript dataprocessing/process_data.r
	# source $(CONDAROOT)/bin/activate && conda activate Rmap-daily-update && Rscript process_radiation_fluxes.r
	source $(CONDAROOT)/bin/activate && conda activate Rmap-daily-update && python dataprocessing/traffic/create_traffic_matrix.py
	source $(CONDAROOT)/bin/activate && conda activate Rmap-daily-update && python dataprocessing/traffic/process_alt_traffic.py
	# source $(CONDAROOT)/bin/activate && conda activate Rmap-daily-update && python dataprocessing/process_regions.py data/uk_forward_commute_flow.csv data/uk_reverse_commute_flow.csv data/areas.csv data/region-groupings.json --threshold=0.8
	source $(CONDAROOT)/bin/activate && conda activate Rmap-daily-update && python dataprocessing/process_regions.py data/traffic_flux_row-normed.csv data/traffic_flux_transpose_row-normed.csv data/areas.csv data/region-groupings.json --threshold=0.8
