//
// Created by ward_s on 5/22/18.
//

#include "../include/spinw.h"
#include "../include/ndmat.h"

magRot spinw::magstr(bool isExact, arma::urowvec newExt, arma::rowvec newOrigin) {
    magRot thisMagRot = magRot();

    matom_struct theseMatoms = matom();
    arma::uword  nAtom = theseMatoms.r.n_cols;
    arma::uword nK = mag_str1.nK;

    double nExt0 = mag_str1.nMagExt;

    // F is stored in 3 * nAtom*nExt0, nK;

    auto extRep = arma::ceil(newExt/nExt0);
    arma::field<arma::cx_vec>M0(nAtom, extRep, nK);

    if (arma::any(newExt != nExt0)){
    } else {
        int ind = 0;

        // Create M0 Matrix
        for (int i = 0; i < nAtom; i++) {
            // We want to replicate along the nExt dimension.
            for (int k = 0; k < extRep; k++) {
                int thisInd = ind;
                for (int j = 0; j < nExt0; j++) {
                    for (int k = 0; k < nK; k++) {
                        M0(i, j, k) = arma::cx_vec(
                                arma::vec(&(mag_str1.F_real[thisInd]), 3),
                                arma::vec(&(mag_str1.F_imag[thisInd]), 3)
                        );
                        thisInd += 3;
                    }
                }
            }
            ind += 3*nK*nExt0;
        }


    }
}
