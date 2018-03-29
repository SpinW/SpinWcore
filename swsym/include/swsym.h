//
// Created by ward_s on 28/02/18.
//

#ifndef SPINW_SWSYM_H
#define SPINW_SWSYM_H

#include <armadillo>
#include <string>

typedef struct Opinfo
{
    bool ismoved;
    bool opmove ;
} Opinfo;


class swsym{
//    arma::mat symOp;
    int max_sym = 300;
    std::string symStr[300];
    std::string symName[300];
    arma::field<arma::cube>symOp;

public:
    explicit swsym(char* dat_dir);
    void interpretSymString(arma::cube &this_cube, std::string inString);
    std::tuple<arma::mat, arma::mat> position(arma::mat &r0, int symN, double tol);
};

#endif //SPINW_SWSYM_H
