//
// Created by ward_s on 28/02/18.
//

#include "../include/swsym.h"
#include "../../include/templateFuncs.h"
#include <iostream>
#include <fstream>
#include <algorithm>


swsym::swsym(char* dat_dirPt){

    symOp = arma::field<arma::cube>(max_sym);

    // What file do we want to load?
    std:: string datDir(dat_dirPt);
    std::string sym_file("symmetry.dat");

    // Create the filename.
    std::string symFile = datDir + sym_file;

    int loop = 0;
    std::string line; //The line taken from the *.txt source
    std::ifstream symFileStream(symFile); //To read from the *.txt File
    std::string temp;
    if (symFileStream.is_open()) {
        while (! symFileStream.eof() ) //Runs while the file is NOT at the end
        {
            if (loop == max_sym){
                std::cout << "Maximum number of symmetry operations loaded (" << max_sym << ")." << std::endl;
                break;
            }
            getline (symFileStream,line); //Gets a single line
            if (line.size() == 0){
                break; // Check if there is a blank line at the end....
            }
            symName[loop] = line.substr(6,11);
            symName[loop].erase(std::remove_if(temp.begin(), temp.end(), ::isspace), temp.end());
            symStr[loop] = line.substr(19,line.size());
            arma::cube tempCube(3,4,30,arma::fill::zeros);
            interpretSymString(tempCube,line.substr(19,line.size()));
            symOp(loop) = tempCube;
            loop++;
        }
        symFileStream.close();
    } else {
        throw std::invalid_argument(std::string("Unable to open dat file"));
//        std::cout << "Unable to open file:" << std::endl << symFile << std::endl;
    }
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
//    this_cube(arma::span::all, arma::span::all, arma::span(nOp+1,this_cube.n_slices-1)).fill(NAN);
}
