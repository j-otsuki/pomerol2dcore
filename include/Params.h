//
// Created by otsuki on 2019/02/11.
//

#ifndef POMEROL2DCORE_PARAMS_H
#define POMEROL2DCORE_PARAMS_H

#include <string>


class Params{
public:
    int n_orb;
    double beta;

    std::string file_h0;
    std::string file_umat;
    std::string file_states;
    std::string file_eigen;
    std::string file_retained;

    // GF
    bool flag_gf;
    int n_w;
    std::string file_gf;

    // Two-particle GF
    bool flag_vx;
    int n_w2f;
    int n_w2b;
    std::string file_vx;

    // methods
    void read(std::string &filename);
    void print();
};

#endif //POMEROL2DCORE_PARAMS_H
