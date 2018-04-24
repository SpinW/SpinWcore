//
// Created by ward_s on 4/24/18.
//

#include "../include/swsym.h"
#include "../../include/templateFuncs.h"

arma::cube swsym::symOperator(int symNumber){
    arma::cube genOP = getSym(symNumber - 1);
    arma::uword nGen  = genOP.n_slices;
    arma::cube symOp;
    symOp.slice(0)= arma::join_rows(arma::eye(3,3), arma::zeros(3,1));

    genOP.each_slice([symOp](const arma::mat& X){
        arma::mat R0 = X(arma::span(),arma::span(0,2));
        arma::mat T0 = X(arma::span(),3);

        arma::mat R = arma::eye(3,3);
        arma::mat T = arma::zeros(3,1);

        //TODO Change P from 1 to the proper value
        int P = 1;

        for (int i = 0; i < P; i++){
            R %= R0;
            T %= R0;
            T += T0;

        }
    });
    return  symOp;
}