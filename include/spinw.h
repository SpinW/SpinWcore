//
// Created by ward_s on 19/02/18.
//

#ifndef UNTITLED2_SPINW_H
#define UNTITLED2_SPINW_H

#include <armadillo>
#include "sw_structs.h"
#include <string>

class spinw {
    lattice lattice1;
    unit_cell unit_cell1;
    twin twin1;
    matrix matrix1;
    single_ion single_ion1;
    coupling coupling1;
    mag_str mag_str1;
    unit unit1;

private:
    matom_struct mag_cache;
    symop_struct sym_cache;

public:

    explicit spinw(lattice latt, unit_cell cell, twin tw, matrix mat, single_ion si, coupling co, mag_str mag, unit un);
    ~spinw();

    arma::mat arma_spinwave(double* qRange, spinwave_opt options);
    std::tuple<arma::cube, arma::mat> arma_twinq(arma::mat q0, arma::cube rotc);
    void intmatrix(int_matrix &this_matrix, bool fitmode, bool plotmode, bool sortDM, bool zeroC, bool extend, bool conjugate);
    std::tuple<arma::mat, arma::urowvec, arma::urowvec> atom();
    matom_struct matom();
    symop_struct symop();
    double* spinwave(double* qRange, spinwave_opt options);
    magRot magstr(bool isExact, arma::urowvec newExt, arma::rowvec newOrigin);
    double* basisvector(bool norm);
};


#endif //UNTITLED2_SPINW_H
