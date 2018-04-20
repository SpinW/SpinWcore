//
// Created by ward_s on 20/02/18.
//

#ifndef SPINW_SW_ADDITIONS_H
#define SPINW_SW_ADDITIONS_H

#include "armadillo"
#include "sw_structs.h"

double* sw_qscan(double* qLim);
arma::mat arma_sw_qscan(double* qLim);
arma::mat arma_basisvector(lattice lattice1, bool norm);
#endif //SPINW_SW_ADDITIONS_H
