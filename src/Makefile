fname = metapop
data_dir = $(repos)/metapop/data

compile:
	clang++ -Wall -Wextra -Wpedantic -o3 -std=c++17 -I${CONDA_PREFIX}/include -lblas -o $(fname) csv.cpp $(fname).cpp

data:
	python make_data.py $(data_dir)

execute:
	./$(fname) 150 $(data_dir)

run:
	make compile && make execute

chart:
	python plot.py $(data_dir)/outputs

clean:
	rm metapop

full:
	make data && make compile && make execute && make chart

