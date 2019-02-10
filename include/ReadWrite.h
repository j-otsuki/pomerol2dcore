//
// Created by otsuki on 2019/02/10.
//

#ifndef POMEROL2DCORE_READWRITE_H
#define POMEROL2DCORE_READWRITE_H

#include <vector>
#include <string>
#include <fstream>


class ReadDataFile
{
private:
    std::ifstream ifs;
    int n_index;
    int n_val;
    std::vector<int> indices;
    std::vector<double> values;
public:
    ReadDataFile(std::string &filename, int n_index, int n_val);
    bool read_line();
    int get_index(int i);
    double get_val(int i);
    void get_indices(std::vector<int> &indices);
    void get_values(std::vector<double> &values);
};


class WriteDataFile
{
private:
    std::ofstream ofs;

public:
    WriteDataFile(std::string &filename);

    template<typename T>
    void write_vector(std::vector<T> &v) {
        for(int i=0; i<v.size(); i++) {
            ofs << v[i] << std::endl;
        }
    }
};

#endif //POMEROL2DCORE_READWRITE_H
