#include "../include/swsym.h"
#include "../../include/templateFuncs.h"

std::tuple<arma::mat, arma::mat> swsym::position(arma::mat &r0, int symN, double tol) {
    arma::cube thisSym = symOp(symN);
    arma::mat r, aIdx, redThisSym, addThisSym;
    arma::uword nAtom = r0.n_cols;

    redThisSym = thisSym(arma::span(),arma::span(0,2),arma::span());
    addThisSym = thisSym(arma::span(),arma::span(3,3),arma::span());
    int perm[] = {1, 3, 2};
    int order[] = {1, 2};
    for (arma::uword i = 0; i < nAtom; i++){
        arma::mat   mTemp = armaModM(mmat(redThisSym, r0(arma::span(), i), order) + addThisSym, 1);
        arma::cube   rTemp(mTemp.n_rows, mTemp.n_cols, 1);
        rTemp.slice(0) = mTemp;
        rTemp  = permute(rTemp, perm);
        if (rTemp.n_elem > 3){
            rTemp = armaModC(rTemp,1);
            arma::find(rTemp > (1 - tol)) -= 1;
            rTemp.slice(0) = uniquetol(rTemp.slice(0),tol);
        }
        r = arma::join_rows(r, rTemp.slice(0));
    }

    return std::make_tuple(r, aIdx);
}