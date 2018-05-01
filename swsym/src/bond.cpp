//
// Created by ward_s on 17/04/18.
//

#include "../include/swsym.h"
#include "../../include/templateFuncs.h"

//TODO test these cases...

std::tuple<arma::mat, arma::urowvec> swsym::bond(arma::mat r,
                                  arma::mat bv,
                                  arma::colvec bond,
                                  int symN,
                                  double toll){

    double tolDist = 1e-5;

    arma::colvec r1 = r.col(bond(3));
    arma::colvec r2 = r.col(bond(4));
    arma::colvec dl = bond(arma::span(0,2));

    arma::cube thisSym = symOperator(symN);
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
        arma::cube r2TempC(r2Temp.memptr(), r2Temp.n_rows, r2Temp.n_cols, 1);
        arma::cube dlTempC(dlTemp.memptr(), dlTemp.n_rows, dlTemp.n_cols, 1);

        r1TempC = permute(r1TempC, perm);
        r2TempC = permute(r2TempC, perm);
        dlTempC = permute(dlTempC, perm);

        r1New = arma::join_rows(r1New,armaModC(r1TempC, 1).slice(0));
        r2New = arma::join_rows(r2New,armaModC(r2TempC, 1).slice(0));

        arma::cube r1TempC2 = arma::floor(r1TempC);
        arma::cube r2TempC2 = arma::floor(r2TempC);

        r1TempC2.elem(arma::find(arma::abs(r1TempC - r1TempC2) > (1 - toll))) += 1;
        r2TempC2.elem(arma::find(arma::abs(r2TempC - r2TempC2) > (1 - toll))) += 1;

        dlTempC += (-r1TempC2 + r2TempC2);

        dlNew = arma::join_rows(dlNew,dlTempC.slice(0));
    }

    std::tuple<arma::urowvec, arma::urowvec> thisAtom1 = isnewUC(r,r1New,tolDist);
    arma::urowvec atom1Select = std::get<1>(thisAtom1);

    if (arma::any(std::get<0>(thisAtom1))){
        throw std::runtime_error("The generated positions for atom1 are wrong!");
    }

    std::tuple<arma::urowvec, arma::urowvec> thisAtom2 = isnewUC(r,r2New,tolDist);
    arma::urowvec atom2Select = std::get<1>(thisAtom2);
    if (arma::any(std::get<0>(thisAtom2))){
        throw std::runtime_error("The generated positions for atom2 are wrong!");
    }

    arma::mat temp = bv * (r.cols(atom2Select) - r.cols(atom1Select) + dlNew);
    temp %= temp;
    arma::rowvec dist = arma::sqrt(arma::sum(temp, 0));
    arma::urowvec rightDist = (dist - dist(0)) < toll;
    if (!(bool)arma::all(rightDist)){
        std::cout << "Symmetry generated couplings are dropped!" << std::endl;
    }
    arma::mat  atom1Select1 = arma::conv_to<arma::mat>::from(atom1Select);
    arma::mat  atom2Select1 = arma::conv_to<arma::mat>::from(atom2Select);

    arma::mat genCP = arma::join_cols(arma::join_cols(dlNew, atom1Select1), atom2Select1);
    genCP = genCP.cols(rightDist % arma::linspace<arma::urowvec>(0,rightDist.n_elem - 1, rightDist.n_elem));

//    arma::uword nC = genCP.n_cols;
    arma::cube theseBonds(genCP.memptr(), genCP.n_rows,genCP.n_cols,1);


    int perm1[] = {2, 3, 1};
    int perm2[] = {3, 2, 1};

    arma::cube C1 = permute(theseBonds, perm1);
    arma::cube C2 = permute(theseBonds, perm2);

    arma::uword bondCols = theseBonds.n_cols;
    arma::cube nc1(5, bondCols, 1);
    nc1.slice(0) = arma::join_cols(-1*genCP.rows(0,2),(genCP.rows(arma::uvec({4, 3}))));
    nc1 = permute(nc1,perm1);

    arma::umat uniqueB(C1.n_rows,C2.n_cols,arma::fill::zeros);
    for (arma::uword s = 0; s < C1.n_slices; s++){
        uniqueB += ((arma::repmat(C1.slice(s),1,C2.n_cols) != arma::repmat(C2.slice(s),C1.n_rows,1))
                && (arma::repmat(nc1.slice(s),1,C2.n_cols) != arma::repmat(C2.slice(s),nc1.n_rows,1)));
    }

    arma::urowvec uniqueB_r = arma::all(arma::trimatu(uniqueB != 0) || arma::trimatl(arma::umat(uniqueB.n_rows,uniqueB.n_cols,arma::fill::ones)),0);

    return std::make_tuple(genCP,uniqueB_r);
}