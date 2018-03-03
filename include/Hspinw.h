//
// Created by ward_s on 19/02/18.
//

#ifndef SPINW_HSPINW_H
#define SPINW_HSPINW_H


#include "sw_structs.h"

struct Hspinw;
typedef struct Hspinw Hspinw;

struct Hsw_sym;
typedef struct Hsw_sym Hsw_sym;

//spinwave_opt spinwave_opt_ini(spinwave_opt opt);

Hspinw* create_sw(lattice latt, unit_cell cell, twin tw, mag_str mag, unit un, char* dat_dir);
//Hspinw* create_empty_sw();
void destroy_sw(Hspinw *sw);
double* sw_basisvector(Hspinw *sw, bool norm);
double* sw_spinwave(Hspinw* sw, double* qRange, spinwave_opt options);

Hsw_sym* loadsym(char* dat_dir);

#endif //SPINW_HSPINW_H
