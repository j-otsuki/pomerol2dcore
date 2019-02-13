//
// Created by otsuki on 2019/02/11.
//

#ifndef POMEROL2DCORE_PARAMS_H
#define POMEROL2DCORE_PARAMS_H

#include <string>


class Params{
public:
    // mandatory parameters
    int n_orb;
    double beta;
    int flag_spin_conserve;

    // optional
    // 0 [default] : (0, up), (1, up), ..., (0, down), (1, down), ...
    // 1           : (0, up), (0, down), (1, up), (1, down), ...
    int index_order;

    // file names
    std::string file_h0;
    std::string file_umat;
    std::string file_states;
    std::string file_eigen;
    std::string file_retained;

    // GF
    int flag_gf;
    int n_w;
    std::string file_gf;

    // two-particle GF
    int flag_vx;
    int n_w2f;
    int n_w2b;
    std::string file_vx;

    // methods
    void read(std::string &filename);
    void print();
};

#endif //POMEROL2DCORE_PARAMS_H
