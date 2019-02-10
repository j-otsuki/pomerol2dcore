//
// Created by otsuki on 2019/02/11.
//

#include "Params.h"
#include <iostream>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>


void Params::read(std::string &filename)
{
    boost::property_tree::ptree pt;
    boost::property_tree::read_ini(filename, pt);

    n_orb = pt.get<int>("n_orb");
    beta = pt.get<double>("beta");

    file_h0       = pt.get<std::string>("file_h0", "h0.in");
    file_umat     = pt.get<std::string>("file_umat", "umat.in");
    file_states   = pt.get<std::string>("file_states", "states.dat");
    file_eigen    = pt.get<std::string>("file_eigen", "eigenvalues.dat");
    file_retained = pt.get<std::string>("file_retained", "retained.dat");

    flag_gf = pt.get<int>("flag_gf", 1);
    n_w = pt.get<int>("n_w", 1024);
    file_gf = pt.get<std::string>("file_gf", "gf.dat");

    flag_vx = pt.get<int>("flag_vx", 0);
    n_w2f = pt.get<int>("n_w2f", 10);
    n_w2b = pt.get<int>("n_w2b", 1);
    file_vx = pt.get<std::string>("file_vx", "vx.dat");
}

void Params::print(){
    std::cout << " n_orb         = " << n_orb << std::endl;
    std::cout << " beta          = " << beta << std::endl;

    std::cout << " file_h0       = " << file_h0 << std::endl;
    std::cout << " file_umat     = " << file_umat << std::endl;
    std::cout << " file_states   = " << file_states << std::endl;
    std::cout << " file_eigen    = " << file_eigen << std::endl;
    std::cout << " file_retained = " << file_retained << std::endl;

    std::cout << " flag_gf       = " << flag_gf << std::endl;
    std::cout << " n_w           = " << n_w << std::endl;
    std::cout << " file_gf       = " << file_gf << std::endl;

    std::cout << " flag_vx       = " << flag_vx << std::endl;
    std::cout << " n_w2f         = " << n_w2f << std::endl;
    std::cout << " n_w2b         = " << n_w2b << std::endl;
    std::cout << " file_vx       = " << file_vx << std::endl;
}
