# metapop
Discrete time SEIR metapopulation modelling for epidemic spread


## usage
All you should need to do is to run
`conda install -c conda-forge openblas`
to install `openblas`. 
(You should already have clang from conda, when you run `which clang` it should give you the local Rmap one.)

Then to make the toy data run
`make data`
to compile the code type
`make`
then to run it 
`make execute`
and to plot the results
`make chart`.
If all of that is too much typing, just do 
`make full`
instead.

If you want to move the file format to something that is the same
as `uk_cases.csv` then you can do
`make coerce`
and then `/data/metapop/outputs/cases.csv` will be the file you should read.
