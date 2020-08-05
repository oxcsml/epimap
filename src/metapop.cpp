#include <istream>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <iterator>
#include <filesystem>
#include <experimental/iterator>
#include <functional>

#include "xtensor/xarray.hpp"
#include "xtensor/xcsv.hpp"
#include "xtensor/xbuilder.hpp"
#include "xtensor/xview.hpp"
#include "xtensor-blas/xlinalg.hpp"

#include "argparse/argparse.hpp"

#include "csv.h"

using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::ifstream;
using std::ofstream;

namespace fs = std::filesystem;

struct cmd_args {
    double a;
    double gamma;
    xt::xarray<double> initial;
    xt::xarray<double> mobility;
    xt::xarray<double> rt;
    fs::path output_folder;
};

cmd_args parse_cmd(int argc, char** argv);
template <typename dtype> 
void write_csv(const char* fname, const vector<xt::xarray<dtype>> &data);
template <typename dtype> xt::xarray<dtype> read_csv(const string fpath);



/*
 * TODO
 * Read Gostic paper
 *  - from reading the paper it seems like what I am doing is reasonable
 *  - could possibly add in some stochasticity
 *  - It's not clear where they get their estimates of `a` and `gamma` from (I can look at bobby paper)
 *  (Maybe I should look at their refs to see if what we are using to compare methods is defensible?)
 *
 * -> Try out YW town and city example
 * -> Find configurations from gostic example 
 *      - (They just had a specific timeseries of R_t)
 * -> Figure out how to tune the R_t and find out what the right parameters for it should be (speak to bobby)
 * -> run it and speak to YW about the results
 * 
 * Notes:
 * - could later make this simulation stochastic as is done here 
 *      https://github.com/cobeylab/Rt_estimation/blob/master/code/simulation.R
 *      (code from gostic et al)
 */

    //"""
    //Gostic plot the product r_0 * s
    //they also show the infections, which is the derivative of my thing
    //"""
    
int main(int argc, char** argv) {

    auto args = parse_cmd(argc, argv);

    auto initial = args.initial;
    const auto beta = args.rt * args.gamma;
    
    // note that these are views, so will modify the array initial
    auto s = xt::row(initial, 0);
    auto e = xt::row(initial, 1);
    auto i = xt::row(initial, 2);
    auto r = xt::row(initial, 3);
    auto n = xt::zeros_like(xt::row(initial, 0));
    n += s + e + i + r;

    vector<xt::xarray<double>> susceptible = {s};
    vector<xt::xarray<double>> exposed = {e};
    vector<xt::xarray<double>> infected = {i};
    vector<xt::xarray<double>> recovered = {r};
    vector<xt::xarray<double>> population = {n};

    auto emigrate = xt::sum(args.mobility, 1);
    for (unsigned int t=0; t < beta.shape()[0]; ++t) {
        auto st=susceptible[t], et=exposed[t], it=infected[t], rt=recovered[t];
        auto sn = st/n;
        auto bt = xt::row(beta, t);
        s = s -  bt * it * sn + xt::linalg::dot(args.mobility, sn) - emigrate * sn; 
        e += bt * it * sn - args.a * et + xt::linalg::dot(args.mobility, et/n) - emigrate * (et/n); 
        i += args.a * et - args.gamma * it + xt::linalg::dot(args.mobility, it/n) - emigrate * (it/n); 
        r += args.gamma * it + xt::linalg::dot(args.mobility, rt/n) - emigrate * (rt/n); 
        n += xt::linalg::dot(args.mobility, n) - emigrate * n;
       
        susceptible.push_back(s);
        exposed.push_back(e);
        infected.push_back(i);
        recovered.push_back(r);
        population.push_back(n);
    } 

    const vector<string> filenames = {
        "susceptible.csv",
        "exposed.csv",
        "infected.csv",
        "recovered.csv",
        "population.csv" 
    };
    vector<vector<xt::xarray<double>>> all_data = {susceptible, exposed, infected, recovered, population};

    for (vector<string>::size_type i=0; i < filenames.size(); ++i) {
        std::filesystem::path fname = filenames[i];
        auto outpath = args.output_folder / fname;
        write_csv((outpath).c_str(), all_data[i]);
    };
    
    return 0;
} 


template <typename dtype> xt::xarray<dtype> read_csv(const string fpath) {
    ifstream file;
    file.open(fpath);
    auto data =  xt::load_csv<dtype>(file);
    return data;
};


template <typename dtype>
void write_csv(const char* fname, const vector<xt::xarray<dtype>> &data){
    //Data should be the rows of the csv file
    //
    //Took out the header row stuff but maybe will end up putting it back in later?!
    //
    CSVWriter writer(fname);
    for (const auto &vec : data) {
        std::vector<dtype> row(vec.begin(), vec.end());
        writer.write_row(row);
    };
} 


cmd_args parse_cmd(int argc, char** argv) {
    auto to_double = [](const string val){return std::stod(val);};

    argparse::ArgumentParser parser("Metapop: deterministic SEIR metapopulation modelling with time varying R.");
    parser.add_argument("--a")
        .required() .help("Parameter `a` such that mean incubation period is 1/a.")
        .action(to_double);
    parser.add_argument("--gamma")
        .required() 
        .help("Parameter `gamma` such that typical length of time for which a person is infectious is 1/gamma")
        .action(to_double);
    parser.add_argument("--init")
        .required()
        .help("Path to csv containing initial values for S, E, I and R");
    parser.add_argument("--mobility")
        .required()
        .help("Path to csv containing mobility patterns");
    parser.add_argument("--rt")
        .required()
        .help("Path to csv containing R_t timeseries by region."
              " The number of rows in this file will determine the"
              " number of timesteps for the simulation.");
    parser.add_argument("-o", "--output")
        .help("Folder in which to save the output files. Defaults to current directory.")
        .default_value(fs::current_path());

    try {
        parser.parse_args(argc, argv);
    }
    catch (const std::runtime_error& err) {
        std::cout << err.what() << std::endl;
        exit(1);
    }

    cmd_args args;

    args.a = parser.get<double>("--a");
    args.gamma = parser.get<double>("--gamma");
    args.output_folder = parser.get<string>("-o");
   
    auto init_pth = parser.get<string>("--init");
    auto mobility_pth = parser.get<string>("--mobility");
    auto rt_pth = parser.get<string>("--rt");
    args.initial = read_csv<double>(init_pth);
    args. mobility =  read_csv<double>(mobility_pth);
    args.rt =  read_csv<double>(rt_pth);

    return args;
}


