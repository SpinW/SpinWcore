//
// Created by ward_s on 4/25/18.
//

#include "../include/swsym.h"
#include "../../include/templateFuncs.h"

int swsym::oporder(arma::mat symOP){
    arma::mat RN = symOP.cols(0, 2);
    arma::colvec TN = arma::round(symOP.col(3) * 12);

    int N = 0;

    while ((arma::norm(RN - arma::eye(3, 3)) > 10 * std::numeric_limits<double>::epsilon() ||
            arma::norm(TN)
            ) && (N < 10)){
        RN *= symOP.cols(0, 2);
        TN = armaModM(arma::round(symOP.cols(0, 2)*TN + symOP.col(3)), 12);
        N++;
    }
    return N;
}