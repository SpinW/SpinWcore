//
// Created by ward_s on 15/03/18.
//

#include "../include/templateFuncs.h"

template<typename T>
arma::Mat<typename T::elem_type> as_Mat(const T& X)
{
    return {X};
}


template<typename T>
T armaMod(const T &X, int n)
{
    const auto& a = as_Mat(X);
    return a - floor(a/n)*n;
}