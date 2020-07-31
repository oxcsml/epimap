#include <istream>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <valarray>
#include <iterator>
#include <experimental/iterator>
#include "csv.h"
#include "mobility.h"
#include "xtensor/xarray.hpp"
#include "xtensor/xcsv.hpp"
#include "xtensor/xbuilder.hpp"

using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::ifstream;
using std::ofstream;

std::tuple<int, string> parse_cmd(int argc, char** argv);
template <typename dtype> vector<vector<dtype>> valarrays_to_vector(const vector<valarray<dtype>> &arr);
template <typename dtype> void write_csv(
        const char* fname, const vector<vector<dtype>> &data, const vector<string> &names);

// I think these are beta, alpha and gamma respectively
const double attack = 1 / 5. ; // force of attack
const double latent = 1 / 5.; // the mean of this variable is 1 / latent period
const double inf_period = 1 / 14.; // length of period infected before recovering or dying

int main(int argc, char** argv) {

    auto [timesteps, data_dir] = parse_cmd(argc, argv);
    string output_folder = data_dir + "/outputs";
    cout << "output folder is: " << output_folder << endl;
    
    xt::xarray<double> s, e, i, r, n;
    // need some way of initialising the values
    s = {980, 20, 9, 80, 880};
    e = {10, 880,  990, 900, 20};
    i = {10, 100,  1, 20, 100};
    r = {0, 0, 0, 0, 0};
    n = s + e + i + r;

    tuple<xt::xarray<double>> susceptible = {s};
    tuple<xt::xarray<double>> exposed = {e};
    tuple<xt::xarray<double>> infected = {i};
    tuple<xt::xarray<double>> recovered = {r};
    tuple<xt::xarray<double>> population = {n};

    ifstream file;
    file.open(data_dir + "/mobility.csv");
    auto mobility =  xt::load_csv<double>(file);

    for (int t=0; t != timesteps; ++t) {
        //xt::xarray<double> sn = s / n;
        s = s - attack * i * s / n;
        e = e + attack * i * s / n - latent * e;
        i = i + latent * e - inf_period * i;
        r = r + inf_period * i;
        //n = n;

        susceptible.push_back(s);
        exposed.push_back(e);
        infected.push_back(i);
        recovered.push_back(r);
        population.push_back(n);
    } 

    const vector<string> region_names = {"0", "1", "2", "3", "4"};
    const vector<string> filenames = {"susceptible.csv", "exposed.csv", "infected.csv", "recovered.csv", "population.csv"};
    vector<vector<valarray<double>>> all_data = {susceptible, exposed, infected, recovered, population};

    for (vector<string>::size_type i=0; i < filenames.size(); ++i) {
        const vector<vector<double>> data = valarrays_to_vector(all_data[i]);
        write_csv((output_folder + filenames[i]).c_str(), data, region_names);
    };

    cout << "R_eff is " << attack / inf_period << endl;
    
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
// --- Turn these populations into things that interact
// --- find some way of testing it?
// --- add in mobility data
// --- add in time varying R_t
// --- send it on to YW
    

template <typename dtype>
void write_csv(const char* fname, const vector<vector<dtype>> &data, const vector<string> &names){
    //* Write data 'matrix' to csv
    //Data should be the rows of the csv file
    //Names are the headers
    CSVWriter writer(fname, names);
    for (const auto &vec : data) {
        writer.write_row(vec);
    };
} 
