fname = metapop
data_dir = $(repos)/metapop/data

compile:
	clang++ -o3 -std=c++17 -I${CONDA_PREFIX}/include -I./include -lblas -o $(fname) $(fname).cpp \
		#-Wall -Wextra -Wpedantic 

data:
	python make_data.py $(data_dir)

execute:
	./$(fname) \
	--a 0.25 \
	--gamma 0.25 \
	--mobility $(data_dir)/basic-mobility.csv \
	--init $(data_dir)/basic-initial-values.csv \
	--rt $(data_dir)/basic-rt.csv \
    --output $(data_dir)/outputs \

run:
	make compile && make execute

chart:
	python plot.py $(data_dir)/outputs

clean:
	rm metapop

full:
	make data && make compile && make execute && make chart

