//
// Created by ward_s on 19/02/18.
//
#include <iostream>
#include "../include/spinw.h"
#include "../include/sw_additions.h"

double* spinw::basisvector(bool norm) {

    arma::mat thisVector = arma_basisvector(this->lattice1, norm);
    auto* basisVector = new double[3*3];
    memcpy(basisVector,thisVector.memptr(),sizeof(double)*3*3);
    return basisVector;
};

double* spinw::spinwave(double* qRange, spinwave_opt options){
  arma::mat result = spinw::arma_spinwave(qRange, options);
    auto* thisResult = new double[result.n_rows*result.n_cols];
//    memcpy(thisResult,result.memptr(),sizeof(double)*result.n_rows*result.n_cols);
    return thisResult;
};