//
// Created by ward_s on 19/02/18.
//

#include <iostream>
#include "spinw.h"
#include "sw_additions.h"

#define PI 3.14159265

spinw::spinw(lattice latt, unit_cell cell, twin tw, mag_str mag, unit un) {
    lattice1 = latt;
    unit_cell1 =cell;
    twin1 = tw;
    mag_str1 = mag;
    unit1 = un;
}

spinw::~spinw(){
    std::cout << "SpinW Object " << this << " Destroyed" << std::endl;
}


arma::mat spinw::arma_spinwave(double* qRange, spinwave_opt options) {

    // For linear scans create the Q line(s)
    arma::mat hkl = arma_sw_qscan(qRange);
    if (options.tid == -2) {
        std::cout << "hkl" << std::endl;
        std::cout << hkl << std::endl;
    }

    // Size of the extended magnetic unit cell
    double nExt_t[3] = {(double) mag_str1.nExt[0], (double) mag_str1.nExt[1], (double) mag_str1.nExt[2]};
    arma::rowvec nExt(&(nExt_t[0]), 3, false, true);
    if (options.tid == -2) {
        std::cout << "nExt" << std::endl;
        std::cout << nExt << std::endl;
    }

    arma::mat k(&(mag_str1.k[0][0]), 3, mag_str1.nMagExt, false, true);
    if (options.tid == -2) {
        std::cout << "k" << std::endl;
        std::cout << k << std::endl;
    }
    arma::mat km = k.t() % nExt;
    if (options.tid == -2) {
        std::cout << "km" << std::endl;
        std::cout << km << std::endl;
    }
    // Check whether the structure is incommensurate
    bool incomm = arma::approx_equal(km, arma::round(km), "absdiff", options.toll);

    // Transform the momentum values to the new lattice coordinate system
    if (options.tid == -2) {
        std::cout << "qmat" << std::endl;
        std::cout << arma::mat(&(unit1.qmat[0][0]), 3, 3, false, true) << std::endl;
    }

    hkl = arma::mat(&(unit1.qmat[0][0]), 3, 3, false, true) * hkl;
    if (options.tid == -2) {
        std::cout << "hkl in lattice" << std::endl;
        std::cout << hkl << std::endl;
    }
    // Calculates momentum transfer in A^-1 units.
    ;
    arma::mat hklA = 2 * PI * (arma::solve(arma_basisvector(false), hkl));
    if (options.tid == -2) {
        std::cout << "hklA" << std::endl;
        std::cout << hklA << std::endl;
    }

    // Check for 2*km
    float new_toll = 2 * options.toll;
//    arma::umat helical = arma::any(
//            arma::sum(arma::square(arma::abs(arma_mod(arma::abs(2 * km) + new_toll, 1) - new_toll)), 1) > new_toll);
//    if (options.tid == -2) {
//        std::cout << "helical?" << std::endl;
//        std::cout << helical << std::endl;
//    }


    // number of Q points
    arma::uword nHkl0 = arma::size(hkl, 1);
    if (options.tid == -2) {
        std::cout << "Number of Q points" << std::endl;
        std::cout << nHkl0 << std::endl;
    }

    // define Q scans for the twins
    int nTwin = twin1.nTwin;
    if (options.notwin) {
        nTwin = 1;
    }

    // if the single twin has no rotation set param.notwin true
    arma::cube rotc(&(twin1.rotc[0][0][0]), 3, 3, twin1.nTwin, false, true);
    arma::mat rotc1 = rotc.slice(0) - arma::eye(3, 3);
    if (options.tid == -2) {
        std::cout << "Twin rotc cut 0" << std::endl;
        std::cout << rotc1 << std::endl;
        std::cout << "Twin rotc cut 0 norm" << std::endl;
        std::cout << arma::norm(arma::vectorise(rotc1)) << std::endl;
    }
    if ((nTwin == 1) & arma::norm(arma::vectorise(rotc1)) == 0) {
        options.notwin = true;
    }

    arma::cube hkl_cube;
    long int nHkl = nHkl0;
    if (!options.notwin) {
        // In the abc coordinate system of the selected twin the scan is
        // rotated opposite direction to rotC.
        std::pair <arma::cube, arma::mat> twinq = spinw::arma_twinq(hkl, rotc);
        nHkl = nHkl0 * nTwin;
        hkl_cube = twinq.first;
    } else {
        hkl_cube.slice(0) = hkl;
        nHkl = nHkl0;
    }
    if (options.tid == -2) {
        std::cout << "HKL cube" << std::endl;
        // This gives a seg fault for some reason :-(
//        std::cout << (*hkl_cube) << std::endl;
    }

    arma::cube hkl0(hkl_cube), hklExt(hkl_cube);
    if (incomm){
        hkl0.each_slice( [](arma::mat& X){ arma::solve(X, bv) * bv;}, true);
    }




//
//    arma::cube F0R(&(mag_str1.F_real[0][0][0]),3,mag_str1.nMagExt,mag_str1.nK,false,true);
//    arma::cube F0I(&(mag_str1.F_imag[0][0][0]),3,mag_str1.nMagExt,mag_str1.nK,false,true);
//    arma::cube F0(F0R,F0I);
//    arma::cube RF0 = arma::sqrt(arma::sum(F0R%F0R,0));
//    arma::cube IF0 = arma::sqrt(arma::sum(F0I%F0I,0));
//
//    arma::cube e3 = F0R/RF0;
//    arma::cube e1 = F0I/IF0;
//    e1 = e1 - ((arma::sum(e1%e3,1)%e3));


    return hkl;
};

std::pair<arma::cube, arma::mat> spinw::arma_twinq(arma::mat q0, arma::cube rotc){

    static arma::mat bv = spinw::arma_basisvector(false);
    static arma::mat this_q0(q0);

    rotc.each_slice( [](arma::mat& X){ arma::solve(X, bv) * bv;}, true);

    arma::cube Qtwin(rotc);
    Qtwin.each_slice( [](arma::mat& X){this_q0 * X;}, true);

    std::pair <arma::cube, arma::mat> return_pair (Qtwin, rotc);
    return return_pair;
}

arma::mat spinw::arma_basisvector(bool norm){

    double alpha = lattice1.angle[0];
    double beta  = lattice1.angle[1];
    double gamma = lattice1.angle[2];

    arma::vec v1 = {1, 0, 0};
    arma::vec v2 = {cos(gamma), sin(gamma), 0};
    arma::vec v3(3);

    v3(0) = cos(beta);
    v3(1) = sin(beta)*(cos(alpha)-cos(beta)*cos(gamma))/(sin(beta)*sin(gamma));
    v3(2) = sqrt(sin(beta)*sin(beta) - v3(1)*v3(1));

    arma::mat thisVector = arma::join_horiz(arma::join_horiz(v1,v2),v3);

    if (!norm) {
        thisVector = thisVector * arma::diagmat(arma::vec(&(lattice1.lat_const[0]),3,false,true));
    };
    return thisVector;
};