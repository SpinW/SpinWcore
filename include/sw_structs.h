//
// Created by ward_s on 19/02/18.
//

#ifndef SPINW_SW_STRUCTS_H
#define SPINW_SW_STRUCTS_H

typedef struct lattice {
    double angle[3];
    double lat_const[3];
    double sym[3*4*100]; // We have up to 100 symmetry operations..
    double origin[3];
    char* label;
    int nSymOp;
} lattice;

typedef struct unit_cell{
    double r[3*100];// We have up to 100 atoms in the unit cell..
    double S[100];
    char** label;
    int color[3*100];
    double ox[100];
    double occ[100];
    double b[2*100];
    double ff[2*11*100];
    int A[100];
    int Z[100];
    double biso[100];
    int nAtom;
} unit_cell;

typedef struct twin {
    double vol[10]; // We have up to 10 twins
    double rotc[3*3*10];
    int nTwin;
} twin;

typedef struct matrix{
    double mat[3*3*100]; // Up to 100 magnetic atoms
    int color[3*100];
    char** label;
    int nMat;
} matrix;

typedef struct single_ion{
    int aniso[100];
    int g[100];
    double field[3];
    double T;
    int nMagAtom;
} single_ion;

typedef struct coupling{
    double dl[3*1000];
    double atom1[1000];
    double atom2[1000];
    int mat_idx[3*1000];
    double idx[1000];
    int type[3*1000];
    int sym[3*1000];
    double rdip;
    int nSym;
    int nBond;
} coupling;

typedef struct mag_str{
    double F_real[3*1000*10];
    double F_imag[3*1000*10];
    double k[3*10];
    int nExt[3];
    int nMagExt;
    int nK;
} mag_str;

typedef struct unit{
    double kB;
    double muB;
    double mu0;
    char** label;
    double nformula;
    double qmat[3*3];

# ifdef __cplusplus
    unit(): kB(1000*arma::datum::k_evk),
            muB(0.057883818066),
            mu0(201.335431),
            nformula(0)
    {}

#endif
} unit;

typedef struct spinwave_opt{
    bool notwin;
    bool sortMode;
    int optmem;
    float toll;
    bool hermit;
    float omega_toll;
    double formfact;
    bool gtensor;
    bool cmplxBase;
    int tid;
    int fid;

    // We can not initialise in C but can in c++
# ifdef __cplusplus
    spinwave_opt(): notwin(false), sortMode(true), optmem(0), toll(1e-4),
                    hermit(true), omega_toll(1e-5), formfact(false),
                    gtensor(false),cmplxBase(false), tid(-1), fid(-1) {}

#endif
} spinwave_opt;

# ifdef __cplusplus
typedef struct matom_struct{
    arma::mat r;
    arma::urowvec idx;
    arma::rowvec S;
} matom_struct;

typedef struct SS{
    arma::mat iso;
    arma::mat ani;
    arma::mat dm;
    arma::mat gen;
    arma::mat bq;
    arma::mat dip;
    arma::mat all;
} SS;

typedef struct SI{
    arma::cube aniso;
    arma::cube g;
    arma::rowvec field;
} SI;

typedef struct int_matrix{
    struct SS SS;
    struct SI SI;
    arma::mat RR;
} int_matrix;

typedef struct JJstruct{
    arma::uvec idx;
    arma::uvec sym;
    arma::cube mat;
    arma::mat type;
    arma::mat iso;
    arma::mat bq;
    arma::mat ani;
    arma::mat dm;
    arma::mat gen;
} JJstruct;

typename struct magRot{
    arma::mat S;
    arma::rowvec k;
    arma::urowvec n;
    arma::urowvec N_ext;
    bool exact;
} magRot;

#endif

#endif //SPINW_SW_STRUCTS_H
