fname = metapop
data_dir = ../../data/metapop
data_name = metapop
a = 0.25

run:
	python $(fname).py \
	--a $(a) \
	--gamma 0.25 \
	--coupling $(data_dir)/$(data_name)-coupling.csv \
	--init $(data_dir)/$(data_name)-initial-values.csv \
	--rt $(data_dir)/$(data_name)-rt.csv \
    --output $(data_dir)/outputs

data:
	python make_data.py $(data_dir)

coerce:
	python coerce.py \
		$(data_dir)/outputs/exposed.csv \
		--output=$(data_dir)/outputs/cases.csv\
		--shift \
		--scale $(a)

csv:
	make run && make coerce

chart:
	python plot.py $(data_dir)/outputs --rt=$(data_dir)/$(data_name)-rt.csv --cases=$(data_dir)/outputs/cases.csv

full:
	make data && make csv && make chart
	
