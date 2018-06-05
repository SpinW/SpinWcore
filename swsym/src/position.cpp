#include "../include/swsym.h"
#include "../../include/templateFuncs.h"

std::tuple<arma::mat, arma::urowvec> swsym::position(arma::mat &r0, int symN, double tol) {
    /*
     * `[r, aIdx] = swsym.position(r0,sym,tol)`
     * generates all symmetry equivalent atomic positions from a given space group and
     * coordinates of the symmetry inequivalent atoms.
     *
     * This is equivalent to calling
     * swsym.position(swsym.generator(SYM),R0)
     */

    arma::cube thisSym = symOp(symN);
    std::tuple<arma::mat, arma::urowvec> returnTP = positionSym(r0, thisSym, tol);

    return returnTP;
}

std::tuple<arma::mat, arma::urowvec, opInfo> swsym::fullPosition(arma::mat &r0, int symN, double tol) {
    /*
     * `[r, aIdx] = swsym.position(r0,sym,tol)`
     * generates all symmetry equivalent atomic positions from a given space group and
     * coordinates of the symmetry inequivalent atoms.
     *
     * This is equivalent to calling
     * swsym.position(swsym.generator(SYM),R0)
     */

    arma::cube thisSym = symOp(symN);
    std::tuple<arma::mat, arma::urowvec, opInfo> returnTP = fullPositionSym(r0, thisSym, tol);

    return returnTP;
}

std::tuple<arma::mat, arma::urowvec> swsym::positionSym(arma::mat &r0, arma::cube &thisSym, double tol) {
    std::tuple<arma::mat, arma::urowvec, opInfo> tempGet = fullPositionSym(r0, thisSym, tol);
    return std::make_tuple(std::get<0>(tempGet), std::get<1>(tempGet));
}


std::tuple<arma::mat, arma::urowvec, opInfo> swsym::fullPositionSym(arma::mat &r0, arma::cube &thisSym, double tol) {
    /*
     * `[r, aIdx] = swsym.position(r0,sym,tol)`
     * generates all symmetry equivalent atomic positions from a given space group and
     * coordinates of the symmetry inequivalent atoms.
     */

    arma::mat r;
    arma::uword nAtom = r0.n_cols;
    arma::urowvec aIdx;

    opInfo thisOpInfo = opInfo();
    thisOpInfo.isMoved = arma::field<arma::urowvec>(nAtom);
    thisOpInfo.opMove = arma::cube(3, 3, 0);

    arma::cube redThisSym = thisSym(arma::span(), arma::span(0, 2), arma::span());
    arma::cube addThisSym = thisSym(arma::span(), arma::span(3, 3), arma::span());

    int order[] = {1, 2};

    for (arma::uword i = 0; i < nAtom; i++) {
        arma::mat thisRTemp;

        for (arma::uword s = 0; s < thisSym.n_slices; s++) {
            arma::mat mTemp = armaModM(mmat(redThisSym.slice(s), r0(arma::span(), i), order) + addThisSym.slice(s), 1);
            thisRTemp = arma::join_rows(thisRTemp, mTemp);
        }
        arma::mat temp = r0.col(i);
        temp = cmodM(temp, tol);
        arma::mat temp2 = cmodM(thisRTemp, tol);
        temp2.each_col([temp](arma::colvec &thisCol){
            thisCol -= temp;
        });
        thisOpInfo.isMoved(i) = arma::sum(arma::pow(temp2, 2),0) > tol*tol;
        arma::urowvec idxF;

        if (thisRTemp.n_elem > 3) {
            thisRTemp = cmodM(thisRTemp, tol);
            std::tuple<arma::mat, arma::urowvec> RV = uniquetolFI(thisRTemp, tol);
            thisRTemp = std::get<0>(RV);
            idxF = std::get<1>(RV);
        } else {
             idxF = arma::urowvec(1, arma::fill::zeros);
        }

        r = arma::join_rows(r, thisRTemp);
        // Notice i+1 so we are in the same formalism as SpinW
        aIdx = arma::join_rows(aIdx, arma::ones<arma::urowvec>(thisRTemp.n_cols) * (i+1));

        thisOpInfo.opMove.insert_slices(thisOpInfo.opMove.n_slices, cubeSlice(redThisSym, idxF, false));
    }

    return std::make_tuple(r, aIdx, thisOpInfo);
}

std::tuple<arma::mat, arma::urowvec> swsym::positionExpand(arma::mat &r0, int symN, double tol) {
    /*
     * `[r, aIdx] = swsym.position(r0,sym,tol)`
     * generates all symmetry equivalent atomic positions from a given space group and
     * coordinates of the symmetry inequivalent atoms.
     *
     * This is equivalent to calling
     * swsym.position(swsym.operator(swsym.generator(SYM)),R0)
     */

    arma::cube thisSym = symOperator(symN);
    std::tuple<arma::mat, arma::urowvec> returnTP = positionSym(r0, thisSym, tol);
    return returnTP;
}

std::tuple<arma::mat, arma::urowvec, opInfo> swsym::fullPositionExpand(arma::mat &r0, int symN, double tol) {
    /*
     * `[r, aIdx] = swsym.position(r0,sym,tol)`
     * generates all symmetry equivalent atomic positions from a given space group and
     * coordinates of the symmetry inequivalent atoms.
     *
     * This is equivalent to calling
     * swsym.position(swsym.operator(swsym.generator(SYM)),R0)
     */

    arma::cube thisSym = symOperator(symN);
    std::tuple<arma::mat, arma::urowvec, opInfo> returnTP = fullPositionSym(r0, thisSym, tol);
    return returnTP;
}