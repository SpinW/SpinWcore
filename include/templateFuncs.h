//
// Created by ward_s on 15/03/18.
//

#ifndef SPINW_TEMPLATEFUNCS_H
#define SPINW_TEMPLATEFUNCS_H

#include <armadillo>
template <typename T> const arma::Mat<T>& as_Mat(const arma::Mat<T>& m) {return m;}
template <typename T> T armaMod(const T &X, int n);
template <typename T> static arma::Cube<T> permute(arma::Cube<T>& cube, const std::tuple<arma::uword,arma::uword,arma::uword>& order);

#endif //SPINW_TEMPLATEFUNCS_H
