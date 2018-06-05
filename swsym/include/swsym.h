//
// Created by ward_s on 28/02/18.
//

#ifndef SPINW_SWSYM_H
#define SPINW_SWSYM_H

#include <armadillo>
#include <string>

typedef struct opInfo {
    arma::field<arma::urowvec> isMoved;
    arma::cube opMove;
} opInfo;


class swsym{

public:
    int totalSyms = 0;
    int max_sym = 300;

protected:
    arma::field<arma::cube>symOp;
    std::string symStr[300];
    std::string symName[300];


public:
    explicit swsym();
    void addSymString(std::string symStr);
    void interpretSymString(arma::cube &this_cube, std::string inString);
    std::tuple<arma::mat, arma::urowvec> position(arma::mat &r0, int symN, double tol);
    std::tuple<arma::mat, arma::urowvec, opInfo> fullPosition(arma::mat &r0, int symN, double tol);
    std::tuple<arma::mat, arma::urowvec> positionExpand(arma::mat &r0, int symN, double tol);
    std::tuple<arma::mat, arma::urowvec, opInfo> fullPositionExpand(arma::mat &r0, int symN, double tol);
    std::tuple<arma::mat, arma::urowvec> bond(arma::mat r, arma::mat bv, arma::colvec bond, int symN, double toll);
    arma::cube getSym(int ind){
        return symOp(ind);
    }

    int searchSym(std::string searchString);
    arma::cube symOperator(int symNumber);
    int oporder(arma::mat symOP);
private:
    std::tuple<arma::mat, arma::urowvec> positionSym(arma::mat &r0, arma::cube &thisSym, double tol);
    std::tuple<arma::mat, arma::urowvec, opInfo> fullPositionSym(arma::mat &r0, arma::cube &thisSym, double tol);

};


#endif //SPINW_SWSYM_H
