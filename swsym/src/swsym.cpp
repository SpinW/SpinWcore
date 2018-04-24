//
// Created by ward_s on 28/02/18.
//

#include "../include/swsym.h"
#include "../../include/templateFuncs.h"
#include "../../include/symmetryDat.h"
#include <iostream>
#include <fstream>
#include <algorithm>


swsym::swsym() {

    ALL_SYM_DAT loadedSyms = ALL_SYM_DAT();
    symOp = arma::field<arma::cube>(max_sym);

    for (int i = 0; i < loadedSyms.symStrLen; i++) {
        if (totalSyms == max_sym) {
            throw std::out_of_range("Maximum number of symmetry operations loaded");
        }
        std::string line = loadedSyms.allStrings[i];
        if (line.size() == 0) {
            break; // Check if there is a blank line at the end....
        }
        symName[totalSyms] = line.substr(6, 11);
        symStr[totalSyms] = line.substr(19, line.size());
        arma::cube tempCube(3, 4, 30, arma::fill::zeros);
        interpretSymString(tempCube, line.substr(19, line.size()));
        symOp(totalSyms) = tempCube;
        totalSyms++;
    }
}

void swsym::addSymString(std::string symStr){
    arma::cube tempCube(3,4,30,arma::fill::zeros);
    interpretSymString(tempCube,symStr);
    char result[100];
    sprintf(result,"Added_%d",totalSyms);
    symName[totalSyms] = std::string(result);
    symOp(totalSyms) = tempCube;
    totalSyms++;
}

int swsym::searchSym(std::string searchString){
    int index = -1;
    for(int x = 0; x < totalSyms; x++){
        if (symName[x].find(searchString, 0) != std::string::npos){
            index = x;
            break;
        }
    }
    return index;
}

void swsym::interpretSymString(arma::cube &this_cube, std::string symStr) {
    int ii = 0;
    int nNew  = 0;
    int nOp   = 0;
    int nSign = 1;

    arma::rowvec vNew(3,arma::fill::zeros);

    while (ii<=symStr.size())
    {
        std::string this_string = symStr.substr(ii,1);

        if (this_string == ",")
        {
            this_cube(arma::span(nNew,nNew), arma::span(0,2), arma::span(nOp,nOp)) = vNew;
            vNew.fill(0);
            nSign = 1 ;
            nNew = (int)(3*((double)nNew/3 - floor((double)nNew/3)) +1);
        }
        else if (this_string  == ";")
        {
            this_cube(arma::span(nNew,nNew), arma::span(0,2), arma::span(nOp,nOp)) = vNew;
            vNew.fill(0);
            nSign = 1 ;
            nNew = 0 ;
            nOp = nOp+1 ;
        }
        else if (this_string  == "x")
        {
            vNew(0) = nSign ;
        }
        else if (this_string  == "y")
        {
            vNew(1) = nSign ;
        }
        else if (this_string == "z")
        {
            vNew(2) = nSign ;
        }
        else if (this_string == "-")
        {
            nSign = -1;
        }
        else if (this_string == "+")
        {
            nSign = 1;
        }
        else if (this_string == "1" || this_string == "2" || this_string == "3")
        {
            this_cube(nNew, 3, nOp) = (double)(symStr[ii]-'0')/(double)(symStr[ii+2]-'0');
            ii+=2;
        }
        ii++;
    }
    this_cube(arma::span(nNew,nNew), arma::span(0,2), arma::span(nOp,nOp)) = vNew ;
    this_cube.shed_slices(nOp+1,this_cube.n_slices-1);
}
