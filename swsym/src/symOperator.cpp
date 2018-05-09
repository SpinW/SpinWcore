//
// Created by ward_s on 4/24/18.
//

#include "../include/swsym.h"
#include "../../include/templateFuncs.h"

inline arma::ucolvec operator|(const arma::colvec& lhs, arma::colvec& rhs){
    // Operate on an element by element basis, and return
    // a urowvec. Decide on something reasonable if the vectors
    // differ in size.

    return arma::any(arma::join_rows(lhs,rhs), 1);
}

arma::cube swsym::symOperator(int symNumber){

    arma::cube genOP = symOp(symNumber);

    genOP.each_slice([](arma::mat &X){
        X.col(3) = arma::round(X.col(3)*12);
    });

    arma::cube symOp(3,4,1);
    symOp.slice(0)= arma::join_rows(arma::eye(3,3), arma::zeros(3,1));

    for (arma::uword s = 0; s < genOP.n_slices; s++){
        arma::mat X = genOP.slice(s);
        arma::mat R0 = X.cols(0,2);
        arma::colvec T0 = X.col(3);

        int P = oporder(arma::join_rows(R0, T0/12));
        arma::mat R = arma::eye(3,3);
        arma::colvec T = arma::zeros(3,1);

        for (int i = 0; i < P; i++){

            R = R0*R;
            T = R0*T + T0;
            arma::uword nSym = symOp.n_slices;

            for (arma::uword j = 0; j < nSym; j++){
                arma::mat temp = symOp.slice(j);
                arma::mat RS = R * temp.cols(0, 2);
                arma::mat TS = armaModM(arma::round(R*(12*temp.col(3)) + T), 12);

                arma::colvec idxR(symOp.n_slices);
                arma::colvec idxT(symOp.n_slices);
                for (arma::uword ss = 0; ss < symOp.n_slices; ss++){
                    idxR(ss) = arma::accu(arma::abs(symOp.slice(ss).cols(0,2) - RS));
                    idxT(ss) = arma::any((12*(symOp.slice(ss).col(3)) - TS));
                }

                if (arma::all(idxR | idxT)){
                    symOp = arma::join_slices(symOp,arma::join_rows(RS, TS/12));
                }
            }
        }
    }
    return  symOp;
}