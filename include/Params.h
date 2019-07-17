//
// pomerol2dcore
//
// Copyright (C) 2019 Junya Otsuki
//

#ifndef POMEROL2DCORE_PARAMS_H
#define POMEROL2DCORE_PARAMS_H

#include <string>


class Params{
public:
    // mandatory parameters
    unsigned short n_orb;
    double beta;
    bool flag_spin_conserve;

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
    std::string file_occup;

    // GF
    bool flag_gf;
    unsigned int n_wf;
    std::string file_gf;

    // susceptibility
    bool flag_suscep;
    int n_wb;
//    std::string file_suscep;
    std::string dir_suscep;

    // two-particle GF
    bool flag_vx;
    int n_w2f;
    int n_w2b;
//    std::string file_vx;
    std::string dir_vx;
    std::string file_freqs;

    // methods
    void read(std::string &filename);
    void print();
};

#endif //POMEROL2DCORE_PARAMS_H
