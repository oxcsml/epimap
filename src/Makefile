fname = metapop
data_dir = $(repos)/metapop/data
data_name = gostic-single

compile:
	clang++ -o3 -std=c++17 -I${CONDA_PREFIX}/include -I./include -lblas -o $(fname) $(fname).cpp \
		#-Wall -Wextra -Wpedantic 

data:
	python make_data.py $(data_dir)

execute:
	./$(fname) \
	--a 0.25 \
	--gamma 0.25 \
	--mobility $(data_dir)/$(data_name)-mobility.csv \
	--init $(data_dir)/$(data_name)-initial-values.csv \
	--rt $(data_dir)/$(data_name)-rt.csv \
    --output $(data_dir)/outputs && \
	python coerce.py $(data_dir)/outputs/infected.csv --output=$(data_dir)/outputs/cases.csv

run:
	make compile && make execute

chart:
	python plot.py $(data_dir)/outputs --rt=$(data_dir)/$(data_name)-rt.csv

clean:
	rm metapop

full:
	make data && make compile && make execute && make chart

