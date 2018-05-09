//
// Created by ward_s on 19/02/18.
//

#include <iostream>
#include "../include/spinw.h"
#include "../include/sw_additions.h"

#define PI 3.14159265

spinw::spinw(lattice latt, unit_cell cell, twin tw, matrix mat, single_ion si, coupling co, mag_str mag, unit un) {
    lattice1 = latt;
    unit_cell1 = cell;
    twin1 = tw;
    matrix1 = mat;
    single_ion1 = si;
    coupling1 = co;
    mag_str1 = mag;
    unit1 = un;
}

spinw::~spinw(){
    std::cout << "SpinW Object " << this << " Destroyed" << std::endl;
}

std::tuple<arma::cube,arma::mat> spinw::arma_twinq(arma::mat q0, arma::cube rotc) {

    std::cout << q0 << std::endl;
    static arma::mat bv = arma_basisvector(this->lattice1, false);
    std::cout << bv << std::endl;
    static arma::mat this_q0(q0);
    std::cout << this_q0 << std::endl;

    rotc.each_slice([](arma::mat &X) { arma::solve(X, bv) * bv; }, true);

    arma::cube Qtwin(rotc);
    Qtwin.each_slice([](arma::mat &X) { this_q0 * X; }, true);

    std::tuple<arma::cube, arma::mat> return_pair(Qtwin, rotc);
    return return_pair;
}
