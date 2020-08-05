#include <fstream>
#include <vector>
#include <experimental/iterator>


class CSVWriter {
    private:
        std::ofstream file;

    public:
        template <class vec> void write_row(const vec &data);
        void write_header(const std::vector<std::string> &names);
        const char delimiter[3] = ",";
        const char newline[2] = "\n";

    CSVWriter(const char* fname) : file(fname) {};
};

template <class vec> void CSVWriter::write_row(const vec &data){
    std::copy(data.begin(), data.end(), std::experimental::make_ostream_joiner(file, delimiter));
    file << newline;
} 

