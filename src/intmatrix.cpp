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
    arma::cube thisMat(&(matrix1.mat[0]),3, 3, matrix1.nMat, false, false);

    arma::uword nMat = (arma::uword)matrix1.nMat;
    thisMat.insert_slices(nMat, arma::zeros<arma::cube>(3, 3, 1));
    thisMat.insert_slices(nMat + 1, arma::cube(3,3,1));
    thisMat.slice(nMat + 1) = 2*arma::eye(3, 3);

    arma::Row<int>thisAniso(&(single_ion1.aniso[0]), single_ion1.nMagAtom);
    arma::Row<int>thisG(&(single_ion1.g[0]), single_ion1.nMagAtom);
    arma::rowvec thisField(&(single_ion1.field[0]), 3);

    int perm[3] = {2, 1, 3};

    // anisotropy matrix
    if(thisAniso.n_cols == nMagAtom){
        arma::urowvec redAniso = arma::conv_to<arma::urowvec>::from(thisAniso); // Comes direct from MATLAB so starts from 1
        redAniso.replace(0, nMat + 1);
//        redAniso.elem(arma::find(redAniso == 0)) += (nMat + 1); // Hence we add one here
        redAniso -= 1;  // Subtract 1 so that we're in sensible indexing.
        if (redAniso.n_elem == 1) {
            this_matrix.SI.aniso = arma::cube(thisMat.n_rows, thisMat.n_cols, 1);
            this_matrix.SI.aniso.slice(0) = thisMat.slice(redAniso(0));
        } else {
            this_matrix.SI.aniso = cubeSlice(thisMat,redAniso);
        }
        this_matrix.SI.aniso += permute(this_matrix.SI.aniso, perm);
        this_matrix.SI.aniso /= 2;
    } else {
        this_matrix.SI.aniso = arma::cube(3, 3, nMagAtom, arma::fill::zeros);
    }

    // g-tensor
    if(thisG.n_cols == nMagAtom){
        arma::urowvec redG = arma::conv_to<arma::urowvec>::from(thisG);
//        redG.elem(arma::find(redG == 0)) += (nMat + 2);
        redG.replace(0, nMat + 2);
        redG -= 1;  // Note the indexing is -1;
        if (redG.n_elem == 1) {
            this_matrix.SI.g = arma::cube(thisMat.n_rows, thisMat.n_cols, 1);
            this_matrix.SI.g.slice(0) = thisMat.slice(redG(0));
        } else {
            this_matrix.SI.g = cubeSlice(thisMat,redG);
        }
        this_matrix.SI.g += permute(this_matrix.SI.g, perm);
        this_matrix.SI.g /= 2;
    } else {
        this_matrix.SI.g = arma::cube(3, 3, nMagAtom);
        this_matrix.SI.g.each_slice([](arma::mat &M){
            M = 2 * arma::mat(M.n_rows, M.n_cols, arma::fill::eye);
        });
    }

    // Bonds
    arma::mat allBonds(&coupling1.dl[0], 3, coupling1.nBond,false,false);

    allBonds.insert_rows(allBonds.n_rows,
                         arma::rowvec(&(coupling1.atom1[0]),coupling1.nBond) - 1);

    allBonds.insert_rows(allBonds.n_rows,
                         arma::rowvec(&(coupling1.atom2[0]),coupling1.nBond) - 1);

    arma::rowvec couplingIDX(&(coupling1.idx[0]),coupling1.nBond);
    allBonds.insert_rows(allBonds.n_rows,couplingIDX - 1);

    this_matrix.SS.all = allBonds;

    arma::uvec lastSym = arma::find(couplingIDX <= (double)coupling1.nSym, 1, "last");


    // extract the assigned bonds

    //TODO maybe this is coupling1.nsym Need to check......

    arma::umat mat_idx = arma::conv_to<arma::umat>::from(
            arma::Mat<int>(&(coupling1.mat_idx[0]),3,coupling1.nBond,false,false));
    arma::umat mat_idxT = mat_idx.t();

    arma::mat mat_type = arma::conv_to<arma::mat>::from(
            arma::Mat<int>(&(coupling1.type[0]),3,coupling1.nBond,false,false));
    arma::mat mat_typeT = mat_type.t();

    arma::umat mat_sym = arma::conv_to<arma::umat>::from(
            arma::Mat<int>(&(coupling1.sym[0]),3,coupling1.nBond,false,false));
    arma::umat mat_symT = mat_sym.t();

    JJstruct JJ = JJstruct();
    JJ.idx = arma::vectorise(mat_idxT.elem(arma::find(mat_idxT != 0))) - 1;
    JJ.sym = arma::vectorise(mat_symT.elem(arma::find(mat_idxT != 0)));

    // keep the column index of each generated bond
    arma::ucolvec colSel  = arma::join_cols(
            arma::find(mat_idx.row(0) != 0),
            arma::join_cols(
                    arma::find(mat_idx.row(1) != 0),
                    arma::find(mat_idx.row(2) != 0)
            )
    );

    this_matrix.SS.all  = this_matrix.SS.all.cols(colSel);
    this_matrix.SS.all.insert_rows(this_matrix.SS.all.n_rows,
                                   arma::conv_to<arma::rowvec>::from(JJ.idx.t()));
    arma::mat toAdd = mat_typeT.elem(arma::find(mat_idx != 0)).t();
    this_matrix.SS.all.insert_rows(this_matrix.SS.all.n_rows, toAdd);

    JJ.mat = arma::cube(thisMat.n_rows, thisMat.n_cols, JJ.idx.n_elem);
    for (arma::uword i = 0; i < JJ.idx.n_elem; i ++){
        JJ.mat.slice(i) = thisMat.slice(JJ.idx(i));
    }
//    std::cout << (this_matrix.SS.all) << std::endl << (JJ.mat) << std::endl;


}