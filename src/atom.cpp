//
// Created by ward_s on 13/03/18.
//

#include "../include/atom.h"
#include "../swsym/include/swsym.h"

std::tuple<arma::mat, arma::urowvec, arma::urowvec> spinw::atom() {
 swsym thisSym = swsym();
 arma::mat thisR(&(unit_cell1.r[0][0]),3,unit_cell1.nAtom);
 std::tuple<arma::mat, arma::urowvec> thesePos = thisSym.position(thisR,lattice1.nSymOp,1E-4);

 arma::urowvec theseInds = std::get<1>(thesePos);

 arma::rowvec S(&(unit_cell1.S[0]),unit_cell1.nAtom);
 arma::urowvec magIDS = S(theseInds) > 0;

 return std::make_tuple(std::get<0>(thesePos),theseInds,magIDS);
 }