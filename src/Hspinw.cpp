//
// Created by ward_s on 19/02/18.
//

#include "../include/spinw.h"
#include "../swsym/include/swsym.h"
#include <iostream>

extern "C" {
    #include "../include/Hspinw.h"
}

Hspinw* create_sw(lattice lat, unit_cell cell, twin tw, mag_str mag, unit un){
    try {
        return reinterpret_cast<Hspinw*>(new spinw(lat, cell, tw, mag, un));
    } catch (const std::exception &e) {
        std::cout << __func__ << "This has all gone wrong" << std::endl;
    }
}

void destroy_sw(Hspinw* sw){
    try{
        delete reinterpret_cast<spinw*>(sw);
    }catch (const std::exception &e) {
        std::cout << __func__ << ": This has all gone wrong" << std::endl;
        std::cout << e.what() << std::endl;
    }
}

double* sw_basisvector(Hspinw* sw, bool norm){
    try{
        return reinterpret_cast<spinw*>(sw)->basisvector(norm);
    }catch (const std::exception &e) {
        std::cout << __func__ << ": This has all gone wrong" << std::endl;
    }
}

double* sw_spinwave(Hspinw* sw, double* qRange, spinwave_opt options){
    try{
        return reinterpret_cast<spinw*>(sw)->spinwave(qRange, options);
    }catch (const std::exception &e) {
        std::cout << __func__ << ": This has all gone wrong" << std::endl;
        std::cout << e.what() << std::endl;
    }
}

Hsw_sym* loadsym(){
    try{
        return reinterpret_cast<Hsw_sym*>(new swsym());
    } catch (const std::exception &e) {
        std::cout << __func__ << ": This has all gone wrong" << std::endl;
        std::cout << e.what() << std::endl;
    }
}