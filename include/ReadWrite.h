//
// Created by otsuki on 2019/02/10.
//

#ifndef POMEROL2DCORE_READWRITE_H
#define POMEROL2DCORE_READWRITE_H

#include <vector>
#include <string>
#include <fstream>
#include <ostream>
#include <complex>


// overload output of complex class
template <typename T>
std::ostream& operator<<(std::ostream& out, const std::complex<T>& c)
{
    out << c.real() << " " << c.imag();
    return out;
}


// A class for reading formatted data in a file
class ReadDataFile
{
private:
    std::ifstream ifs;
    unsigned int n_index;
    unsigned int n_val;
    std::vector<int> indices;
    std::vector<double> values;
public:
    ReadDataFile(std::string &filename, unsigned int n_index, unsigned int n_val);
    bool read_line();
    int get_index(int i);
    double get_val(int i);
    void get_indices(std::vector<int> &indices);
    void get_values(std::vector<double> &values);
};


// A class for writing data into a file
class WriteDataFile
{
private:
    std::ofstream ofs;

public:
    explicit WriteDataFile(std::string &filename);

    template<typename T>
    void write_vector(const std::vector<T> &v) {
        for(int i=0; i<v.size(); i++) {
            ofs << v[i] << std::endl;
        }
    }
    void write_str(const std::string &str) {
        ofs << str << std::endl;
    }
};

#endif //POMEROL2DCORE_READWRITE_H
