//
// Created by ward_s on 20/02/18.
//

#include "../include/sw_additions.h"

//template<typename T> T arma_mod(T a, int n)
//{
//    return a - n*arma::floor(a/n);
//}


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

