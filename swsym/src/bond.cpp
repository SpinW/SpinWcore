//
// Created by ward_s on 17/04/18.
//

#include "../include/swsym.h"
#include "../../include/templateFuncs.h"

std::tuple<arma::mat, arma::umat> swsym::bond(arma::mat r,
                                  const arma::mat &bv,
                                  arma::rowvec bond,
                                  int symN,
                                  double toll){

    double tolDist = 1e-5;

    arma::colvec r1 = r(arma::span(),bond(3));
    arma::colvec r2 = r(arma::span(),bond(4));
    arma::rowvec dl = bond(arma::span(0,2));

    arma::cube thisSym = symOp(symN);
    arma::cube redThisSym = thisSym(arma::span(), arma::span(0, 2), arma::span());
    arma::cube addThisSym = thisSym(arma::span(), arma::span(3, 3), arma::span());

    int perm[] = {1, 3, 2};
    int order[] = {1, 2};

    arma::mat r1New, r2New, dlNew;
    for (arma::uword s = 0; s < thisSym.n_slices; s++) {

        arma::mat r1Temp = mmat(redThisSym.slice(s), r1, order) + addThisSym.slice(s);
        arma::mat r2Temp = mmat(redThisSym.slice(s), r2, order) + addThisSym.slice(s);
        arma::mat dlTemp = mmat(redThisSym.slice(s), dl, order) + addThisSym.slice(s);

        arma::cube r1TempC(r1Temp.memptr(), r1Temp.n_rows, r1Temp.n_cols, 1);
        arma::cube r2TempC(r1Temp.memptr(), r2Temp.n_rows, r2Temp.n_cols, 1);
        arma::cube dlTempC(dlTemp.memptr(), dlTemp.n_rows, dlTemp.n_cols, 1);
//
//        r1TempC.slice(0) = r1Temp;
//        r2TempC.slice(0) = r2Temp;
//        dlTempC.slice(0) = dlTemp;

        r1TempC = permute(r1TempC, perm);
        r2TempC = permute(r2TempC, perm);
        dlTempC = permute(dlTempC, perm);

        r1TempC = armaModC(r1TempC, 1); auto r1TempC2 = r1TempC;
        r2TempC = armaModC(r2TempC, 1); auto r2TempC2 = r2TempC;

        r1TempC2.elem(arma::find(r1TempC2 > (1 - toll))) -= 1;
        r2TempC2.elem(arma::find(r2TempC2 > (1 - toll))) -= 1;

        dlTempC += (-r1TempC2 + r2TempC2);

        dlTempC = armaModC(dlTempC, 1);

        r1New = arma::join_rows(r1New,r1TempC.slice(0));
        r2New = arma::join_rows(r2New,r2TempC.slice(0));
        dlNew = arma::join_rows(dlNew,dlTempC.slice(0));
    }
    std::tuple<arma::uvec, arma::uvec> atom1 = isnewUC(r,r1New,tolDist);
    arma::uvec atom1Select = std::get<1>(atom1);

    if (arma::any(std::get<0>(atom1))){
        throw std::runtime_error("The generated positions for atom1 are wrong!");
    }

    std::tuple<arma::uvec, arma::uvec> atom2 = isnewUC(r,r2New,tolDist);
    arma::uvec atom2Select = std::get<1>(atom2);
    if (arma::any(std::get<0>(atom2))){
        throw std::runtime_error("The generated positions for atom2 are wrong!");
    }

    arma::mat temp = (r.cols(atom2Select) - r.cols(atom1Select) + dlNew);
    temp = temp % temp;
    arma::vec dist = arma::sqrt(arma::sum(temp, 0));
    arma::uvec rightDist = (dist - dist(0)) < toll;
    if (!(bool)arma::all(rightDist)){
        std::cout << "Symmetry generated couplings are dropped!" << std::endl;
    }
    arma::mat  atom1Select1 = arma::conv_to<arma::mat>::from(atom1Select);
    arma::mat  atom2Select1 = arma::conv_to<arma::mat>::from(atom2Select);

    arma::mat genCP = arma::join_rows(arma::join_rows(dlNew, atom1Select1), atom2Select1);
    genCP = genCP.rows(rightDist);

    arma::uword nC = genCP.n_cols;
    arma::cube theseBonds(genCP.memptr(), genCP.n_rows,genCP.n_cols,1);
//    theseBonds.slice(0) = genCP;

    int perm1[] = {2, 3, 1};
    int perm2[] = {3, 2, 1};
    arma::cube C1 = permute(theseBonds, perm1);
    arma::cube C2 = permute(theseBonds, perm2);

    arma::uword bondCols = theseBonds.n_cols;
    arma::cube nc1(bondCols,5,1);
    nc1.slice(0) = arma::join_cols(genCP.rows(0,2),genCP.rows(3,4));
    nc1 = permute(nc1,perm1);

    arma::umat uniqueB(C1.n_rows,C1.n_cols,arma::fill::zeros);
    for (arma::uword s = 0; s < C1.n_slices; s++){
        uniqueB += (C1.slice(s) != C2.slice(s));
    }
    uniqueB = arma::trimatu(uniqueB) || arma::trimatu(arma::umat(uniqueB.n_rows,uniqueB.n_cols,arma::fill::ones));

    return std::make_tuple(genCP,uniqueB);
}