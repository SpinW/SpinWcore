//
// Created by ward_s on 13/03/18.
//

#include "../include/atom.h"
#include "../swsym/include/swsym.h"

std::tuple<arma::mat, arma::urowvec, arma::urowvec> spinw::atom() {
    // Return r, idx, mag. Note that mag is logical not an index of magnetic.

    swsym thisSym = swsym();
    arma::mat thisR(&(unit_cell1.r[0]),3, unit_cell1.nAtom,false,true);
    std::tuple<arma::mat, arma::urowvec> thesePos = thisSym.positionExpand(thisR, lattice1.nSymOp, 1E-4);

    arma::urowvec theseInds = std::get<1>(thesePos);

    arma::rowvec S(&(unit_cell1.S[0]), unit_cell1.nAtom);
    arma::uvec magIDS = S.elem(theseInds - 1) > 0; // -1 as arrays start from 0;

    return std::make_tuple(std::get<0>(thesePos), theseInds, magIDS.t());
}

matom_struct spinw::matom() {

    if (mag_cache.idx.is_empty()) {

        arma::rowvec thisS(&(unit_cell1.S[0]), unit_cell1.nAtom);

        // Return r, idx, S.
        std::tuple<arma::mat, arma::urowvec, arma::urowvec> theseAtoms = atom();

        arma::mat thisR = std::get<0>(theseAtoms); // R
        arma::urowvec thisIdx = std::get<1>(theseAtoms); //IDX Atom Type number
        arma::urowvec magIDS = std::get<2>(theseAtoms); //Is Spin (logical)

        if (magIDS.n_elem > 1) {
            mag_cache.r = thisR.cols(arma::find(magIDS));
            mag_cache.idx = thisIdx.elem(arma::find(std::get<2>(theseAtoms))).t();
            mag_cache.S = thisS.cols(mag_cache.idx - 1);
        } else {
            mag_cache.r = thisR;
            mag_cache.idx = thisIdx;
            mag_cache.S = thisS;
        }
    }
    return mag_cache;
}