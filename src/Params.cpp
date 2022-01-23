//
// pomerol2dcore
//
// Copyright (C) 2019 Junya Otsuki
//

#include "Params.h"
#include <iostream>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>


void Params::read(std::string &filename)
{
    boost::property_tree::ptree pt;
    boost::property_tree::read_ini(filename, pt);

    // mandatory parameters
    n_orb = pt.get<unsigned short>("n_orb");
    beta = pt.get<double>("beta");
    flag_spin_conserve = pt.get<bool>("flag_spin_conserve");

    // optional
    n_bath = pt.get<unsigned int>("n_bath", 0);
    index_order = pt.get<int>("index_order", 0);
    tol_real = pt.get<double>("tol_real", 1e-12);

    // file names
    file_h0       = pt.get<std::string>("file_h0", "h0.in");
    file_umat     = pt.get<std::string>("file_umat", "umat.in");
    file_states   = pt.get<std::string>("file_states", "states.dat");
    file_eigen    = pt.get<std::string>("file_eigen", "eigenvalues.dat");
    file_retained = pt.get<std::string>("file_retained", "retained.dat");
    file_occup    = pt.get<std::string>("file_occup", "occup.dat");

    // GF
    flag_gf = pt.get<bool>("flag_gf", true);
    n_wf = pt.get<unsigned int>("n_wf", 1024);
    file_gf = pt.get<std::string>("file_gf", "gf.dat");

    // susceptibility
    flag_suscep = pt.get<bool>("flag_suscep", false);
    n_wb = pt.get<unsigned int>("n_wb", 100);
    dir_suscep = pt.get<std::string>("dir_suscep", "susceptibility");

    // two-particle GF
    flag_vx = pt.get<bool>("flag_vx", false);
    n_w2f = pt.get<int>("n_w2f", 10);
    n_w2b = pt.get<int>("n_w2b", 1);
    dir_vx = pt.get<std::string>("dir_vx", "two_particle");
    file_freqs = pt.get<std::string>("file_freqs", "");
}

void Params::print()
{
    std::cout << " -------------------------------- Mandatory" << std::endl;
    std::cout << " n_orb               = " << n_orb << std::endl;
    std::cout << " beta                = " << beta << std::endl;
    std::cout << " flag_spin_conserve  = " << flag_spin_conserve << std::endl;

    std::cout << " -------------------------------- Optional" << std::endl;
    std::cout << " n_bath              = " << n_bath << std::endl;
    std::cout << " index_order         = " << index_order << std::endl;
    std::cout << " tol_real            = " << tol_real << std::endl;

    std::cout << " -------------------------------- File names" << std::endl;
    std::cout << " file_h0             = " << file_h0 << std::endl;
    std::cout << " file_umat           = " << file_umat << std::endl;
    std::cout << " file_states         = " << file_states << std::endl;
    std::cout << " file_eigen          = " << file_eigen << std::endl;
    std::cout << " file_retained       = " << file_retained << std::endl;
    std::cout << " file_occup          = " << file_occup << std::endl;

    std::cout << " -------------------------------- GF" << std::endl;
    std::cout << " flag_gf             = " << flag_gf << std::endl;
    std::cout << " n_wf                = " << n_wf << std::endl;
    std::cout << " file_gf             = " << file_gf << std::endl;

    std::cout << " -------------------------------- Susceptibility" << std::endl;
    std::cout << " flag_suscep         = " << flag_suscep << std::endl;
    std::cout << " n_wb                = " << n_wb << std::endl;
    std::cout << " dir_suscep          = " << dir_suscep << std::endl;

    std::cout << " -------------------------------- Two-particle GF" << std::endl;
    std::cout << " flag_vx             = " << flag_vx << std::endl;
    std::cout << " n_w2f               = " << n_w2f << std::endl;
    std::cout << " n_w2b               = " << n_w2b << std::endl;
    std::cout << " dir_vx              = " << dir_vx << std::endl;
    std::cout << " file_freqs          = " << file_freqs << std::endl;
}
