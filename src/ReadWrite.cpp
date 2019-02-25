//
// pomerol2dcore
//
// Copyright (C) 2019 Junya Otsuki
//

#include "ReadWrite.h"
#include <iostream>
#include <sstream>
#include <cassert>
#include <iomanip>


ReadDataFile::ReadDataFile(std::string &filename, unsigned int n_index, unsigned int n_val) : n_index(n_index), n_val(n_val) {
    ifs.open(filename);
    if( ifs.fail() ){
        std::cerr << "Failed in opening the file" << std::endl;
        exit(2);
    }
    indices.resize(n_index);
    values.resize(n_val);
}

// return true if succeeded
bool ReadDataFile::read_line() {
    bool status = false;
    std::string line;
    if( std::getline(ifs, line) ){
        // split by white space
        std::stringstream ss(line);
        std::string word;
        std::vector<std::string> words;
        while( std::getline(ss, word, ' ') ){
            words.push_back(word);
        }
        assert( words.size() == n_index + n_val );
        // set indices
        for(int i=0; i<n_index; i++){
            indices[i] = std::stoi(words[i]);
        }
        // set values
        for(int i=0; i<n_val; i++){
            values[i] = std::stod(words[i+n_index]);
        }
        status = true;
    }
    return status;
}

int ReadDataFile::get_index(int i) {
    assert(i < n_index);
    return indices[i];
}

double ReadDataFile::get_val(int i) {
    assert(i < n_val);
    return values[i];
}

void ReadDataFile::get_indices(std::vector<int> &indices) {
    indices = this->indices;
}

void ReadDataFile::get_values(std::vector<double> &values) {
    values = this->values;
}


WriteDataFile::WriteDataFile(std::string &filename){
    ofs.open(filename);
    if( ofs.fail() ){
        std::cerr << "Failed in opening the file" << std::endl;
        exit(2);
    }
    ofs << std::scientific << std::setprecision(15);
}
