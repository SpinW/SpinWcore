//
// Created by ward_s on 01/05/18.
//

#include "../include/intmatrix.h"
#include "../include/templateFuncs.h"


void spinw::intmatrix(struct init_matrix &this_matrix, bool fitmode, bool plotmode, bool sortDM, bool zeroC, bool extend, bool conjugate){

    arma::rowvec nExt0 = {(double)mag_str1.nExt[0],(double)mag_str1.nExt[1],(double)mag_str1.nExt[2]};

    if (arma::sum(nExt0) == 3){
        extend = false;
    }

    matom_struct mAtom = mag_cache;
    arma::uword nMagAtom = mAtom.r.n_cols;
    arma::cube thisMat(&(matrix1.mat[0][0][0]), 3, 3, matrix1.nMat);
    arma::uword nMat = matrix1.nMat;
    thisMat.insert_slices(thisMat.n_slices,arma::zeros<arma::mat>(3, 3));
    thisMat.insert_slices(thisMat.n_slices,arma::eye<arma::mat>(3, 3));

    arma::Row thisAniso<int>(&(single_ion1.aniso[0]),single_ion1.nMagAtom);

    arma::Row thisG<int>(&(single_ion1.g[0]),single_ion1.nMagAtom);
    arma::rowvec thisField(&(single_ion1.field[0]),3);

    int perm[3] = {2, 1, 3};

    // anisotropy matrix
    if(thisAniso.n_cols == nMagAtom){
        thisAniso.elem(arma::find(thisAniso == 0)) = nMat + 1;
        this_matrix.SI.aniso = thisMat.slices(thisAniso - 1); // Note the indexing is -1;
        this_matrix.SI.aniso += permute(this_matrix.SI.aniso, perm);
        this_matrix.SI.aniso /= 2;
    } else {
        this_matrix.SI.aniso = arma::cube(3, 3, nMagAtom, arma::fill::zeros);
    }

    // g-tensor
    if(thisG.n_cols == nMagAtom){
        thisG.elem(arma::find(thisG == 0)) = nMat + 2;
        this_matrix.SI.g = thisMat.slices(thisG - 1); // Note the indexing is -1;
        this_matrix.SI.g += permute(this_matrix.SI.g, perm);
        this_matrix.SI.g /= 2;
    } else {
        this_matrix.SI.g = arma::cube(3, 3, nMagAtom, arma::fill::eye);
        this_matrix.SI.g *= 2;
    }

    // Bonds
    arma::mat allBonds (&(coupling1.dl[0][0]), 3, coupling1.nBond);
    allBonds.insert_rows(allBonds.n_rows, arma::rowvec(&(coupling1.atom1[0]),coupling1.nBond) - 1);
    allBonds.insert_rows(allBonds.n_rows, arma::rowvec(&(coupling1.atom2[0]),coupling1.nBond) - 1);
    allBonds.insert_rows(allBonds.n_rows, arma::rowvec(&(coupling1.idx[0]),coupling1.nBond) - 1);
}