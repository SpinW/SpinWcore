#include "../include/swsym.h"
#include "../../include/templateFuncs.h"

std::tuple<arma::mat, arma::urowvec> swsym::position(arma::mat &r0, int symN, double tol) {
    /*
     * `[r, aIdx] = swsym.position(r0,sym,tol)`
     * generates all symmetry equivalent atomic positions from a given space group and
     * coordinates of the symmetry inequivalent atoms.
     */

    arma::cube thisSym = symOp(symN);
    arma::mat r;
    arma::uword nAtom = r0.n_cols;
    arma::urowvec aIdx;

    arma::cube redThisSym = thisSym(arma::span(), arma::span(0, 2), arma::span());
    arma::cube addThisSym = thisSym(arma::span(), arma::span(3, 3), arma::span());

    int perm[] = {1, 3, 2};
    int order[] = {1, 2};

    for (arma::uword i = 0; i < nAtom; i++) {
        for (arma::uword s = 0; s < thisSym.n_slices; s++) {
            arma::mat mTemp = armaModM(mmat(redThisSym.slice(s), r0(arma::span(), i), order) + addThisSym.slice(s), 1);
            arma::cube rTemp(mTemp.memptr(), mTemp.n_rows, mTemp.n_cols, 1);
            rTemp = permute(rTemp, perm);
            if (rTemp.n_elem > 3) {
                rTemp = armaModC(rTemp, 1);
                rTemp.elem(arma::find(rTemp > (1 - tol))) -= 1;
                rTemp.slice(0) = uniquetol(rTemp.slice(0), tol);
            }
            r = arma::join_rows(r, rTemp.slice(0));
            aIdx = arma::join_rows(aIdx, arma::ones<arma::urowvec>(rTemp.n_cols) * i);
        }
    }

    return std::make_tuple(r, aIdx);
}