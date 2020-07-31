fname = metapop

compile:
	clang++ -Wall -Wextra -Wpedantic -o3 -std=c++17 -I${CONDA_PREFIX}/include -o $(fname) csv.cpp $(fname).cpp

execute:
	./$(fname) 150 $(repos)/metapop/data

run:
	make compile && make execute

chart:
	python plot.py $(repos)/metapop/outputs/

clean:
	rm metapop

data:
	python make_data.py ${repos}/metapop/data
