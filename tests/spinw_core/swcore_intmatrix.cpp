//
// Created by ward_s on 09/05/18.
//

#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "../../src/templateFuncs.tpp"
#include "../../include/spinw.h"

using namespace std;
#define PI 3.14159265

TEST(SwCoreIntMatrix, M1){

    unit_cell thisUnitCell = unit_cell();
    lattice thisLattice = lattice();
    mag_str thisMagStr = mag_str();
    coupling thisCoupling = coupling();
    single_ion thisIon = single_ion();
    matrix thisMat = matrix();

    // Make the Unit Cell
    thisUnitCell.r[0] = 0;
    thisUnitCell.r[1] = 0;
    thisUnitCell.r[2] = 0;
    thisUnitCell.nAtom = 1;
    thisUnitCell.S[0] = 1;
    thisUnitCell.b[0] = 1;
    thisUnitCell.b[1] = 1;
    int ff[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1};
    for (int i = 0; i < 22; i++){
        thisUnitCell.ff[i] = ff[i];
    }
    thisUnitCell.A[0] = -1;
    thisUnitCell.Z[0] = 113;
    thisUnitCell.biso[0] = 0;
    thisUnitCell.ox[0] = 0;
    thisUnitCell.occ[0] = 1;

    // Make the lattice
    thisLattice.angle[0] = PI/2;
    thisLattice.angle[1] = PI/2;
    thisLattice.angle[2] = 2*PI/3;
    thisLattice.lat_const[0] = 3.0;
    thisLattice.lat_const[1] = 3.0;
    thisLattice.lat_const[2] = 9.0;
    thisLattice.origin[0] = 0;
    thisLattice.origin[1] = 0;
    thisLattice.origin[2] = 0;
    thisLattice.nSymOp = 1;
    int eye[] = {1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0};
    for (int i = 0; i < 12; i++){
        thisLattice.sym[i] = eye[i];
    }

    // Make the mag_str
    thisMagStr.nExt[0] = 0;
    thisMagStr.nExt[1] = 0;
    thisMagStr.nExt[2] = 0;
    thisMagStr.k[0] = 1.0/3;
    thisMagStr.k[1] = 1.0/3;
    thisMagStr.k[2] = 0;
    thisMagStr.F_real[0] = 1;
    thisMagStr.F_real[1] = 0;
    thisMagStr.F_real[2] = 0;
    thisMagStr.F_imag[0] = 0;
    thisMagStr.F_imag[1] = 1;
    thisMagStr.F_imag[2] = 0;
    thisMagStr.nK = 1;
    thisMagStr.nMagExt = 1;

    int dl[] = {0, 1, 0, 1, 0, 0, 1, 1, 0, 1,
                -1, 0, 1, 2, 0, 2, 1, 0, 0, 2,
                0, 2, 0, 0, 2, 2, 0, 1, -2, 0,
                2, -1, 0, 1, 3, 0, 2, 3, 0, 3,
                1, 0, 3, 2, 0, 3, 0, 0, 0, 3,
                0, 3, 3, 0, 0, 0, 1, 1, 0, -1,
                1, 1, -1, 0, -1, 1, 1, 0, 1,
                0, 1, 1, 1, 1, 1};

    thisCoupling.rdip = 0;
    thisCoupling.nSym = 0;
    int idxC[] = {1, 1, 1, 2, 2, 2, 3, 3,
                  3, 4, 4, 4, 4, 4, 4, 5,
                  5, 5, 5, 6, 6, 6, 6, 6, 6};
    for (int i = 0; i < 25; i++){
        thisCoupling.atom1[i] = 1;
        thisCoupling.atom2[i] = 1;
        thisCoupling.type[i*3 +0] = 0;
        thisCoupling.type[i*3 +1] = 0;
        thisCoupling.type[i*3 +2] = 0;

        thisCoupling.dl[i*3 + 0] = dl[i*3 + 0];
        thisCoupling.dl[i*3 + 1] = dl[i*3 + 1];
        thisCoupling.dl[i*3 + 2] = dl[i*3 + 2];

        thisCoupling.idx[i] = idxC[i];

        if (i < 3){
            thisCoupling.sym[i*3 + 0] = 1;
            thisCoupling.sym[i*3 + 1] = 0;
            thisCoupling.sym[i*3 + 2] = 0;
            thisCoupling.mat_idx[i*3 +0] = 1;
            thisCoupling.mat_idx[i*3 +1] = 0;
            thisCoupling.mat_idx[i*3 +2] = 0;
        } else{
            thisCoupling.sym[i*3 + 0] = 0;
            thisCoupling.sym[i*3 + 1] = 0;
            thisCoupling.sym[i*3 + 2] = 0;
            thisCoupling.mat_idx[i*3 +0] = 0;
            thisCoupling.mat_idx[i*3 +1] = 0;
            thisCoupling.mat_idx[i*3 +2] = 0;
        }
    }
    thisCoupling.nBond = 25;

    //Make single ion
    thisIon.aniso[0] = 0;
    thisIon.g[0] = 0;
    thisIon.field[0] = 0;
    thisIon.field[1] = 0;
    thisIon.field[2] = 0;
    thisIon.T = 0;
    thisIon.nMagAtom = 1;

    //Make Matrix
    thisMat.nMat = 1;
    thisMat.mat[0] = 1; thisMat.mat[1] = 0; thisMat.mat[2] = 0;
    thisMat.mat[3] = 0; thisMat.mat[4] = 1; thisMat.mat[5] = 0;
    thisMat.mat[6] = 0; thisMat.mat[7] = 0; thisMat.mat[8] = 1;

    spinw s = spinw(thisLattice, thisUnitCell, twin(), thisMat, thisIon, thisCoupling, thisMagStr, unit());

    int_matrix intMatrixResult = int_matrix();
    s.intmatrix(intMatrixResult, true, false, false, false, false, true);

}