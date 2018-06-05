//
// Created by ward_s on 14/05/18.
//

#include "../include/symop.h"
#include "../include/sw_additions.h"

symop_struct spinw::symop() {

    coupling thisCoupling = coupling1;
    arma::mat bondList(6, thisCoupling.nBond);
    bondList.rows(0, 3) = arma::mat(&(thisCoupling.dl[0]), 3, thisCoupling.nBond);
    bondList.row(3) = arma::rowvec(&(thisCoupling.atom1[0]), thisCoupling.nBond);
    bondList.row(4) = arma::rowvec(&(thisCoupling.atom2[0]), thisCoupling.nBond);
    bondList.row(5) = arma::rowvec(&(thisCoupling.idx[0]), thisCoupling.nBond);

    arma::urowvec lastSym = arma::find(bondList.row(5) <= thisCoupling.nSym, 1, "last");

    if (sym_cache.bond.is_empty()){

        if (thisCoupling.nSym > 0){
            arma::mat BV = arma_basisvector(lattice1, false);
            swsym thisSym = swsym();
//            std::tuple<arma::mat, arma::urowvec, opInfo> thisPos = thisSym.fullPosition(arma::mat &r0, int lattice1., double tol);
        }


    }

    return sym_cache;
}