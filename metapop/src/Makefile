fname = metapop
data_dir = ../../data/metapop
data_name = gostic-single

python:
	python $(fname).py \
	--a 0.25 \
	--gamma 0.25 \
	--mobility $(data_dir)/$(data_name)-mobility.csv \
	--init $(data_dir)/$(data_name)-initial-values.csv \
	--rt $(data_dir)/$(data_name)-rt.csv \
    --output $(data_dir)/outputs

data:
	python make_data.py $(data_dir)

coerce:
	python coerce.py $(data_dir)/outputs/infected.csv --output=$(data_dir)/outputs/cases.csv

csv:
	make python && make coerce

chart:
	python plot.py $(data_dir)/outputs --rt=$(data_dir)/$(data_name)-rt.csv

full:
	make data && make python && make chart
	
