//
// Created by ward_s on 28/02/18.
//

#ifndef SPINW_SWSYM_H
#define SPINW_SWSYM_H

//#include "../include/position.h"
#include <armadillo>

#include <string>


class swsym{
//    arma::mat symOp;
    int max_sym = 300;
    std::string symStr[300];
    std::string symName[300];
    arma::field<arma::cube>symOp;

public:
    explicit swsym(char* dat_dir);
    void interpretSymString(arma::cube &this_cube, std::string inString);
//    void position(TYPE symOp, TYPE r0, TYPE fid, TYPE tol, TYPE& r, TYPE& aIdx, _Opinfo& opInfo) ;

};


#endif //SPINW_SWSYM_H
