#include <fstream>
#include <iostream>
#include <vector>
#include <experimental/iterator>


using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::ofstream;

void print_vector(vector<int> &vec);
void write_csv(const char* fname, const vector<vector<int>> data, const vector<string> names);

class CSVWriter; 

class CSVWriter {
    private:
        std::ofstream file;

    public:
        template <class vec> void write_row(const vec &data);
        void write_header(const vector<std::string> &names);
        void flush();
        void close();
        const char delimiter[3] = ", ";
        const char newline[2] = "\n";

    CSVWriter(const char* fname) : file(fname) {};
};

void CSVWriter::write_header(const vector<std::string> &names) {
    file << "index" << delimiter;
    write_row(names);
} 

template <class vec> void CSVWriter::write_row(const vec &data){
    std::copy(data.begin(), data.end(), std::experimental::make_ostream_joiner(file, delimiter));
    file << newline;
} 

const float attack = 1 / 5. ; // force of attack
const float latent = 1 / 5.; // the mean of this variable is 1 / latent period
const float inf_period = 1 / 14.; // length of period infected before recovering or dying

int main() {
    const int timesteps = 100;
    const int nppl = 1000;

    vector<int> infected  = {10}; 
    vector<int> susceptible = {980};
    vector<int> exposed = {10}; 
    vector<int> recovered  = {0}; 
   
    for (int t=0; t != timesteps; t++) {
        susceptible.push_back(susceptible[t] - attack * infected[t] * susceptible[t] / nppl);
        exposed.push_back(exposed[t] + attack * infected[t] * susceptible[t] / nppl - latent * exposed[t]);
        infected.push_back(infected[t] + latent * exposed[t] - inf_period * infected[t]);
        recovered.push_back(recovered[t] + inf_period * infected[t]);
    } 
  

    const vector<vector<int>> data = {susceptible, exposed, infected, recovered};
    const vector<string> names = {"susceptible", "exposed", "infected", "recovered"};
    
    const vector<int> test_data = {0, 1, 2, 3, 4, 5};
    CSVWriter writer("test.csv");
    writer.write_header(names);
    writer.write_row(test_data);

    write_csv("output.csv", data, names);
    
    //cout << "R_eff is " << attack / inf_period << endl;
    //print_vector(susceptible);
    //print_vector(infected);
    
    return 0;
} 

void print_vector(vector<int> &vec) {
    for (const int& i : vec) 
        std::cout << i << ' ';
    std::cout << '\n';
} 


void write_csv(const char* fname, const vector<vector<int>> data, const vector<string> names){
/** Write data to csv
   Will write only the number of rows corresponding to the length of the shortest
   vector given
*/

    vector<vector<int>::size_type> lengths;
    for (const auto &vec: data) {
        lengths.push_back(vec.size());
    } 

    vector<int>::size_type  min_length = *std::min_element(lengths.begin(), lengths.end());

    CSVWriter writer(fname);
    writer.write_header(names);

    for (unsigned int t=0; t < min_length; ++t) {
        vector<int> row;
        row.push_back(t);
        for (const auto &vec: data) {
            row.push_back(vec[t]);
        };
        writer.write_row(row);
    };
} 
