//
// Created by ward_s on 01/05/18.
//

#include "../include/intmatrix.h"
#include "../include/templateFuncs.h"


void spinw::intmatrix(struct init_matrix &this_matrix, bool fitmode, bool plotmode, bool sortDM, bool zeroC, bool extend, bool conjugate){

    arma::rowvec nExt0 = {(double)mag_str1.nExt[0],
                          (double)mag_str1.nExt[1],
                          (double)mag_str1.nExt[2]};

    if (arma::sum(nExt0) == 3){
        extend = false;
    }

    matom_struct mAtom = matom();
    arma::uword nMagAtom = mAtom.r.n_cols;
    arma::cube thisMat(&(matrix1.mat[0]),3, 3, matrix1.nMat + 2, false, false);

    arma::uword nMat = (arma::uword)matrix1.nMat;
    thisMat.slice(nMat) = arma::zeros<arma::mat>(3, 3);
    thisMat.slice(nMat + 1) = arma::eye<arma::mat>(3, 3);

    arma::Row<int>thisAniso(&(single_ion1.aniso[0]), single_ion1.nMagAtom);
    arma::Row<int>thisG(&(single_ion1.g[0]), single_ion1.nMagAtom);
    arma::rowvec thisField(&(single_ion1.field[0]), 3);

    int perm[3] = {2, 1, 3};

    // anisotropy matrix
    if(thisAniso.n_cols == nMagAtom){
        arma::urowvec redAniso = arma::conv_to<arma::urowvec>::from(thisAniso);
        redAniso.elem(arma::find(redAniso == 0)) += (nMat + 1);
        redAniso -= 1;  // Note the indexing is -1;
        this_matrix.SI.aniso = cubeSlice(thisMat,redAniso);
        this_matrix.SI.aniso += permute(this_matrix.SI.aniso, perm);
        this_matrix.SI.aniso /= 2;
    } else {
        this_matrix.SI.aniso = arma::cube(3, 3, nMagAtom, arma::fill::zeros);
    }

    // g-tensor
    if(thisG.n_cols == nMagAtom){
        arma::urowvec redG = arma::conv_to<arma::urowvec>::from(thisG);
        redG.elem(arma::find(redG == 0)) += (nMat + 2);
        redG -= 1;  // Note the indexing is -1;
        this_matrix.SI.g = cubeSlice(thisMat,redG);
        this_matrix.SI.g += permute(this_matrix.SI.g, perm);
        this_matrix.SI.g /= 2;
    } else {
        this_matrix.SI.g = arma::cube(3, 3, nMagAtom, arma::fill::eye);
        this_matrix.SI.g *= 2;
    }

    // Bonds
//    arma::mat allBonds = arma::join_cols(
//            arma::rowvec(&(coupling1.dl[0][0]),coupling1.nBond,false,false),
//            arma::join_cols(
//                    arma::rowvec(&(coupling1.dl[1][0]),coupling1.nBond,false,false),
//                    arma::rowvec(&(coupling1.dl[2][0]),coupling1.nBond,false,false)
//            )
//    );
    arma::mat allBonds(&coupling1.dl[0],coupling1.nBond,false,false);

    allBonds.insert_rows(allBonds.n_rows,
                         arma::rowvec(&(coupling1.atom1[0]),coupling1.nBond) - 1);
    allBonds.insert_rows(allBonds.n_rows,
                         arma::rowvec(&(coupling1.atom2[0]),coupling1.nBond) - 1);

    arma::rowvec couplingIDX(&(coupling1.idx[0]),coupling1.nBond);
    allBonds.insert_rows(allBonds.n_rows,couplingIDX - 1);

    this_matrix.SS.all = allBonds;


    arma::urowvec lastSym = arma::find(couplingIDX <= (double)coupling1.nSym, 1, "last");


    // extract the assigned bonds

    //TODO maybe this is coupling1.nsym Need to check......

//    arma::mat mat_idx = arma::conv_to<arma::mat>::from(
//            (arma::Mat<int>)arma::join_cols(
//                    arma::Row<int>(&(coupling1.mat_idx[0][0]),coupling1.nBond,false,false),
//            arma::join_cols(
//                    arma::Row<int>(&(coupling1.mat_idx[1][0]),coupling1.nBond,false,false),
//                    arma::Row<int>(&(coupling1.mat_idx[2][0]),coupling1.nBond,false,false)
//            )
//    ));
    arma::mat mat_idx = arma::conv_to<arma::mat>::from(arma::Mat<int>(&(coupling1.mat_idx[0]),coupling1.nBond,false,false));
    arma::mat mat_idxT = mat_idx.t();

//    arma::mat mat_type = arma::conv_to<arma::mat>::from(
//            (arma::Mat<int>)arma::join_cols(
//            arma::Row<int>(&(coupling1.type[0][0]),coupling1.nBond,false,false),
//            arma::join_cols(
//                    arma::Row<int>(&(coupling1.type[1][0]),coupling1.nBond,false,false),
//                    arma::Row<int>(&(coupling1.type[2][0]),coupling1.nBond,false,false)
//            )
//    ));
    arma::mat mat_type = arma::conv_to<arma::mat>::from(arma::Mat<int>(&(coupling1.type[0]),coupling1.nBond,false,false));
    arma::mat mat_typeT = mat_type.t();

//    arma::mat mat_sym = arma::conv_to<arma::mat>::from(
//            (arma::Mat<int>)arma::join_cols(
//                    arma::Row<int>(&(coupling1.sym[0][0]),coupling1.nBond,false,false),
//                    arma::join_cols(
//                            arma::Row<int>(&(coupling1.sym[1][0]),coupling1.nBond,false,false),
//                            arma::Row<int>(&(coupling1.sym[2][0]),coupling1.nBond,false,false)
//                    )
//            ));
    arma::mat mat_sym = arma::conv_to<arma::mat>::from(arma::Mat<int>(&(coupling1.sym[0]),coupling1.nBond,false,false));
    arma::mat mat_symT = mat_sym.t();

    JJstruct JJ = JJstruct();
    JJ.idx = mat_idxT.elem(arma::find(mat_idxT != 0));
    JJ.sym = mat_symT.elem(arma::find(mat_idxT != 0));

    // keep the column index of each generated bond
//    arma::colvec colSel  = arma::join_cols(
//            arma::find(mat_idx.row(0) != 0),
//            arma::join_cols(
//                    arma::find(mat_idx.row(1) != 0),
//                    arma::find(mat_idx.row(2) != 0)
//            )
//    );
//    SS.all  = SS.all.each_row([colSel](arma::rowvec &R){
//        R = R.elem(colSel);
//    });
}