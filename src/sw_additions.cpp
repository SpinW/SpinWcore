//
// Created by ward_s on 20/02/18.
//

#include "../include/sw_additions.h"

double *sw_qscan(double *qLim) {

    size_t array_offset = 2;

    size_t n = qLim[0];
    size_t step = qLim[1];

    auto nq = (int) qLim[n * step + array_offset];

    arma::mat A;
    if (n > 0) {
        A = arma_sw_qscan(qLim);
    } else {
        A = arma::mat(&(qLim[2]), step, nq, false, true);
    }
    auto *test = new double[step * nq * n];
    memcpy(test, A.memptr(), sizeof(double) * step * nq * n);
    return test;
};

arma::mat arma_sw_qscan(double *qLim) {

    size_t array_offset = 2;
    size_t n = qLim[0];
    size_t step = qLim[1];

    auto nq = (size_t) qLim[n * step + array_offset];
    arma::mat A;
    if (n > 0) {
        A = arma::mat(step, nq * (n - 1));
        for (size_t i = 0; i < step; i++) {
            for (size_t j = 0; j < (n - 1); j++) {
                A(i, arma::span(j * nq, (j + 1) * nq - 1)) =
                        arma::linspace<arma::rowvec>(qLim[i + j * n + array_offset],
                                                     qLim[i + j * n + step + array_offset],
                                                     nq);
            }
        }
    } else {
        A = arma::mat(&(qLim[2]), step, nq, false, true);
    }
    return A;
};

arma::mat arma_basisvector(lattice lattice1, bool norm){

    double alpha = lattice1.angle[0];
    double beta  = lattice1.angle[1];
    double gamma = lattice1.angle[2];

//    arma::vec v1 = {1, 0, 0};
//    arma::vec v2 = {cos(gamma), sin(gamma), 0};
//    arma::vec v3(3);
//
//    v3(0) = cos(beta);
//    v3(1) = sin(beta)*(cos(alpha)-cos(beta)*cos(gamma))/(sin(beta)*sin(gamma));
//    v3(2) = sqrt(sin(beta)*sin(beta) - v3(1)*v3(1));

    arma::mat thisVector = {
            {1, cos(gamma), cos(beta)},
            {0, sin(gamma), sin(beta)*(cos(alpha)-cos(beta)*cos(gamma))/(sin(beta)*sin(gamma))},
            {0, 0,          sqrt(sin(beta)*sin(beta) - (sin(beta)*(cos(alpha)-cos(beta)*cos(gamma))/(sin(beta)*sin(gamma)) * sin(beta)*(cos(alpha)-cos(beta)*cos(gamma))/(sin(beta)*sin(gamma))))}
    };
//    arma::mat thisVector = arma::join_horiz(arma::join_horiz(v1,v2),v3);

    if (!norm) {
        thisVector = thisVector * arma::diagmat(arma::vec(&(lattice1.lat_const[0]),3,false,true));
    };
    return thisVector;
};
