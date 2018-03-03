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
    std::string dat_dir;

public:

    explicit spinw(lattice latt, unit_cell cell, twin tw, mag_str mag, unit un, char *dat_folder);
    ~spinw();

    arma::mat arma_spinwave(double* qRange, spinwave_opt options);
    arma::mat arma_basisvector(bool norm);
    std::pair<arma::cube, arma::mat> arma_twinq(arma::mat q0, arma::cube rotc);
    void initmatrix(init_matrix *this_matrix, bool fitmode, bool plotmode, bool sortDM, bool zeroC, bool extend, bool conjugate);
    void matom(matom_struct *this_matom);

    double* spinwave(double* qRange, spinwave_opt options);
    double* basisvector(bool norm);
};


#endif //UNTITLED2_SPINW_H
