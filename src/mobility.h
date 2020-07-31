#include <valarray>
#include <vector>

// This is a bit of a weird implementation of a matrix,
// would be good to change it to something more sensible later on

template <typename dtype>
dtype inner(const std::valarray<dtype> &u, const std::valarray<dtype> &v) {
    assert(u.size() == v.size());
    return (u * v).sum();
}

template <typename dtype>
class MobilityMatrix {
    private:
        std::vector<std::valarray<dtype>> rows;

    public:
        MobilityMatrix(std::vector<std::valarray<dtype>> &data) : rows(data) {};
        std::valarray<dtype> vecmulr(const std::valarray<dtype> &vec);

};

template <typename dtype>
std::valarray<dtype> MobilityMatrix<dtype>::vecmulr(const std::valarray<dtype> &vec) {
    std::vector<dtype> temp; 
    for (const auto &row: rows) {
        temp.push_back(inner(row, vec));
    }
    std::valarray<dtype> output(temp.data(), temp.size());
    return output;
}

