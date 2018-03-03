//
// Created by ward_s on 28/02/18.
//

#include "../include/swsym.h"
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
    std::cout << symFile << std::endl;


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
            getline (symFileStream,line); //Gets a single line from example.txt
            if (line.size() == 0){
                break; // Check if there is a blank line at the end....
            }
            temp = line.substr(6,11);
            temp.erase(std::remove_if(temp.begin(), temp.end(), ::isspace), temp.end());
            symName[loop] = temp;
            symStr[loop] = line.substr(19,line.size());
            loop++; //Does an increment to the variable 'loop'
        }
        symFileStream.close(); //Closes the file
    } else {
        std::cout << "Unable to open file:" << std::endl << symFile << std::endl;
    }
//    for(loop2 = 0; loop2 < (loop-1); loop2++){ //For loop make to cout the lines stored
//        std::cout << loop2 << std::endl;
//        std::cout << array[loop2] << std::endl;
//    }
}

//arma::cube swsym::interpretSymString(std::string symStr) {
//    int ii = 0;
//    int nNew  = 0;
//    int nOp   = 0;
//    int nSign = 1;
//    arma::vec vNew(3,arma::fill::zeros);
//    while (ii<=symStr.size())
//    {
//        if (symStr.substr(ii,1).compare(","))
//        {
//            symOp(nNew, m2cpp::fspan(1, 1, 3), nOp) = vNew ;
//            vNew = vNew*0 ;
//            nSign = 1 ;
//            nNew = mod(nNew, 3)+1 ;
//        }
//        else if (symStr(ii)==";")
//        {
//            symOp(nNew, m2cpp::fspan(1, 1, 3), nOp) = vNew ;
//            vNew = vNew*0 ;
//            nSign = 1 ;
//            nNew = 1 ;
//            nOp = nOp+1 ;
//        }
//        else if (symStr(ii)=="x")
//        {
//            vNew(1) = nSign ;
//        }
//        else if (symStr(ii)=="y")
//        {
//            vNew(2) = nSign ;
//        }
//        else if (symStr(ii)=="z")
//        {
//            vNew(3) = nSign ;
//        }
//        else if (symStr(ii)=="-")
//        {
//            nSign = -1 ;
//        }
//        else if (symStr(ii)=="+")
//        {
//            nSign = 1 ;
//        }
//        else if ((symStr(ii)=="1") || (symStr(ii)=="2") || (symStr(ii)=="3"))
//        {
//            symOp(nNew, 4, nOp) = (symStr(ii)-"0")/(symStr(ii+2)-"0") ;
//            ii+=2 ;
//        }
//        ii++ ;
//    }
//    symOp(nNew, m2cpp::fspan(1, 1, 3), nOp) = vNew ;
//    symOp = symOp(m2cpp::span<uvec>(0, symOp.n_rows-1), m2cpp::span<uvec>(0, symOp.n_cols-1), m2cpp::fspan(1, 1, nOp)) ;
//
//}
