//
// Created by ward_s on 4/24/18.
//

#include "../include/swsym.h"
#include "../../include/templateFuncs.h"

arma::cube swsym::symOperator(int symNumber){


    arma::cube genOP = getSym(symNumber);

    genOP.each_slice([](arma::mat&X){
        X.col(3) = X.col(3)*12;
        X.col(3) = arma::round(X.col(3));
    });

    arma::cube symOp(3,4,1);
    symOp.slice(0)= arma::join_rows(arma::eye(3,3), arma::zeros(3,1));

    for (arma::uword s = 0; s < genOP.n_slices; s++){
        arma::mat X = genOP.slice(s);
        arma::mat R0 = X.cols(0,2);
        arma::colvec T0 = X.col(3);

        arma::mat R = arma::eye(3,3);
        arma::colvec T = arma::zeros(3,1);

        int P = oporder(arma::join_rows(R0, T0));

        for (int i = 0; i < P; i++){
            R = R0*R;
            T = R0*T + T0;

            arma::uword nSym = symOp.n_slices;
            for (arma::uword j = 0; j < nSym; j++){
                arma::mat temp = symOp.slice(j);
                arma::mat RS = R * temp.cols(0, 2);
                arma::mat TS = armaModM(arma::round(R*temp.col(3) + T), 12);

                arma::cube newCube = symOp;
                newCube.each_slice([RS, TS](arma::mat& Y) {
                    Y(0, 0) = arma::accu(arma::abs(Y.cols(0,2) - RS));
                    Y(0, 1) = arma::any(arma::vectorise(Y.col(3) - TS) > 0);
                });

                arma::colvec idxR = newCube.subcube(0, 0, 0, 0, 0, newCube.n_slices-1);
                arma::colvec idxT = newCube.subcube(0, 1, 0, 0, 1, newCube.n_slices-1);

                bool add = true;
                for (int a = 0; a < idxR.n_elem; a++){
                    if ((idxR(a) == 0) && (idxT(a) == 0)){
                        add = false;
                        break;
                    }
                }

                if (add){
                    symOp = arma::join_slices(symOp,arma::join_rows(RS,TS/12));
                }
            }
        }
    }
    return  symOp;
}