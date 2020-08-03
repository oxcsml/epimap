#include <istream>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <iterator>
#include <filesystem>
#include <experimental/iterator>
#include "csv.h"

#include "xtensor/xarray.hpp"
#include "xtensor/xcsv.hpp"
#include "xtensor/xbuilder.hpp"
#include "xtensor/xview.hpp"
#include "xtensor-blas/xlinalg.hpp"

using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::ifstream;
using std::ofstream;

namespace fs = std::filesystem;

std::tuple<int, string> parse_cmd(int argc, char** argv);
template <typename dtype> 
void write_csv(const char* fname, const vector<xt::xarray<dtype>> &data, const vector<string> &names);

// I think these are beta, alpha and gamma respectively
const double beta = 1/5. ; // force of attack, typical time between contacts is 1/beta
const double a = 1 / 5.; // the average incubation period is 1/a
const double inf_period = 1 / 14.; // 1 / length of period infected before recovering or dying (this is gamma on wiki)

template <typename dtype>
xt::xarray<dtype> load_csv(const string fpath) {
    ifstream file;
    file.open(fpath);
    auto data =  xt::load_csv<dtype>(file);
    return data;
};


int main(int argc, char** argv) {
    auto [timesteps, data_dir] = parse_cmd(argc, argv);
    fs::path data_folder = data_dir;
    fs::path output_folder = data_dir / std::filesystem::path("outputs");
    
    auto mobility = load_csv<double>(data_folder / fs::path("mobility.csv"));
    auto initial = load_csv<double>(data_folder / fs::path("initial-values-seed.csv"));

    cout << "Mobility patterns are:\n" << mobility  << endl;
    cout << "Initial values are:\n" << initial << endl;
    cout << "a is " << 1 / a << endl;
    cout << "Infectious period is " << 1 / inf_period << endl;
    cout << "R_eff is " << beta / inf_period << endl;

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

    auto emigrate = xt::sum(mobility, 0);
    for (int t=0; t != timesteps; ++t) {
        auto st=susceptible[t], et=exposed[t], it=infected[t], rt=recovered[t];
        auto sn = st/n;
        s -= beta * it * sn + xt::linalg::dot(mobility, sn) - emigrate * sn; 
        e += beta * it * sn - a * et + xt::linalg::dot(mobility, et/n) - emigrate * (et/n); 
        i += a * et - inf_period * it + xt::linalg::dot(mobility, it/n) - emigrate * (it/n); 
        r += inf_period * it + xt::linalg::dot(mobility, rt/n) - emigrate * (rt/n); 
        n += xt::linalg::dot(mobility, n) - emigrate * n;

        susceptible.push_back(s);
        exposed.push_back(e);
        infected.push_back(i);
        recovered.push_back(r);
        population.push_back(n);
    } 

    const vector<string> region_names = {"0", "1", "2", "3", "4"};
    const vector<string> filenames = {
        "susceptible.csv",
		"exposed.csv"    ,
		"infected.csv"   ,
		"recovered.csv"  ,
		"population.csv" 
    };
    vector<vector<xt::xarray<double>>> all_data = {susceptible, exposed, infected, recovered, population};

    for (vector<string>::size_type i=0; i < filenames.size(); ++i) {
        std::filesystem::path fname = filenames[i];
        auto outpath = output_folder / fname;
        write_csv((outpath).c_str(), all_data[i], region_names);
    };
    
    return 0;
} 

std::tuple<int, string> parse_cmd(int argc, char** argv) {
    int timesteps;
    string data_dir;
    
    if (argc < 3) {
        throw std::domain_error(
                "Not enough arguments given. Must be given integer number of timesteps and path to output folder."
                );
    } else {
       timesteps = std::stoi(argv[1]);
       data_dir = argv[2];
    } 

    std::tuple<int, string> tup(timesteps, data_dir);
    return tup;
} 

/*
 * Integrate metapopulation modelling and put outputs to csv (maybe can use matrices for this?)
 * Read R_t in from a csv too (make python stuff to look at this)
 * Make some sort of visualisation so it looks like it's working
 * Figure out how to tune the R_t and find out what the right parameters for it should be (speak to bobby)
 * Speak to YW about the results
 *
 */
// Later, use the seir initial values to calculate nppl
// for now ca just put itin here since eventually will habe multiple populations anyway
// can just remove the nppl argument because it's irrelevant!
//
// load beta, alpha, gamma from file
// 
// --- Turn these populations into things that interact
// --- find some way of testing it?
// --- add in mobility data
// --- add in time varying R_t
// --- send it on to YW
    

template <typename dtype>
void write_csv(const char* fname, const vector<xt::xarray<dtype>> &data, const vector<string> &names){
    //Data should be the rows of the csv file
    //Names are the headers
    CSVWriter writer(fname, names);
    for (const auto &vec : data) {
        std::vector<dtype> row(vec.begin(), vec.end());
        writer.write_row(row);
    };
} 
