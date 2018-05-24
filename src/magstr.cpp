//
// Created by ward_s on 5/22/18.
//

#include "../include/spinw.h"

magRot spinw::magstr(bool isExact, arma::urowvec newExt, arma::rowvec newOrigin) {
    magRot thisMagRot = magRot;

    matom_struct theseMatoms = matom();
    arma::uword  nAtom = theseMatoms.r.n_cols;
    arma::uword nK = mag_str1.nK;

    double nExt0 = mag_str1.nMagExt;

    arma::field<arma::cube>
    if (arma::any(newExt != nExt0)){

    }
}
