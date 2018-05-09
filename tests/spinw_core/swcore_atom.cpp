//
// Created by ward_s on 03/05/18.
//

#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "../../src/templateFuncs.tpp"
#include "../../include/spinw.h"

using namespace std;

typedef struct atom_struct{
    arma::mat r;
    arma::urowvec idx;
    arma::urowvec mag;
}atom_struct;

class SwCoreAtomResults {
public:
    matom_struct mag_results;
    atom_struct atom_results;

    int symIDX = -1;
    int noAtoms = -1;

    SwCoreAtomResults(int inSym, int inAtom) {
        symIDX = inSym;
        noAtoms = inAtom;
        switch (inSym) {
            case 0:
                switch (noAtoms){
                    case 1:
                        atom_results.r = arma::mat(3, 1, arma::fill::zeros);
                        atom_results.idx = arma::urowvec(1, arma::fill::ones);
                        atom_results.mag = arma::urowvec(1, arma::fill::ones);
                        mag_results.r = arma::mat(3, 1, arma::fill::zeros);
                        mag_results.idx = arma::urowvec(1, arma::fill::ones);
                        mag_results.S = arma::rowvec(1, arma::fill::ones);
                        break;
                    case 2:
                        atom_results.r = arma::mat({{0, 0.1},{0, 0.1},{0, 0}});
                        atom_results.idx = arma::urowvec({1, 2});
                        atom_results.mag = arma::urowvec(2, arma::fill::ones);
                        mag_results.r = arma::mat({{0, 0.1},{0, 0.1},{0, 0}});
                        mag_results.idx = arma::urowvec({1, 2});
                        mag_results.S = arma::rowvec({1, 0.5});
                        break;
                    case 3:
                        atom_results.r = arma::mat({{0, 0.1},{0, 0.1},{0, 0}});
                        atom_results.idx = arma::urowvec({1, 2});
                        atom_results.mag = arma::urowvec({0, 1});
                        mag_results.r = arma::mat({0.1, 0.1, 0}).t();
                        mag_results.idx = arma::urowvec({2});
                        mag_results.S = arma::rowvec({0.5});
                        break;
                }
                break;
            case 9:
                switch (noAtoms) {
                    case 1:
                        atom_results.r = arma::mat(3, 1, arma::fill::zeros);
                        atom_results.idx = arma::urowvec(1, arma::fill::ones);
                        atom_results.mag = arma::urowvec(1, arma::fill::ones);
                        mag_results.r = arma::mat(3, 1, arma::fill::zeros);
                        mag_results.idx = arma::urowvec(1, arma::fill::ones);
                        mag_results.S = arma::rowvec(1, arma::fill::ones);
                        break;
                    case 2:
                        atom_results.r = arma::mat({{0, 0.1, 0.9, 0.9, 0.1}, {0, 0.1, 0.1, 0.9, 0.9}, {0, 0, 0, 0, 0}});
                        atom_results.idx = arma::urowvec({1, 2, 2, 2, 2});
                        atom_results.mag = arma::urowvec(5, arma::fill::ones);
                        mag_results.r = arma::mat({{0, 0.1, 0.9, 0.9, 0.1}, {0, 0.1, 0.1, 0.9, 0.9}, {0, 0, 0, 0, 0}});
                        mag_results.idx = arma::urowvec({1, 2, 2, 2, 2});
                        mag_results.S = arma::rowvec({1, 0.5, 0.5, 0.5, 0.5});
                        break;
                    case 3:
                        atom_results.r = arma::mat({{0, 0.1, 0.9, 0.9, 0.1}, {0, 0.1, 0.1, 0.9, 0.9}, {0, 0, 0, 0, 0}});
                        atom_results.idx = arma::urowvec({1, 2, 2, 2, 2});
                        atom_results.mag = arma::urowvec({0, 1, 1, 1, 1});
                        mag_results.r = arma::mat({{0.1, 0.9, 0.9, 0.1}, {0.1, 0.1, 0.9, 0.9}, {0, 0, 0, 0}});
                        mag_results.idx = arma::urowvec({2, 2, 2, 2});
                        mag_results.S = 0.5*arma::rowvec(4, arma::fill::ones);
                        break;
                }
                break;
            case 19:
                switch (noAtoms) {
                    case 1:
                        atom_results.r = arma::mat({{0, 0.5, 0,   0.5},
                                                    {0, 0.5, 0,   0.5},
                                                    {0, 0,   0.5, 0.5}});
                        atom_results.idx = arma::urowvec(4, arma::fill::ones);
                        atom_results.mag = arma::urowvec(4, arma::fill::ones);
                        mag_results.r = arma::mat({{0, 0.5, 0,   0.5},
                                                   {0, 0.5, 0,   0.5},
                                                   {0, 0,   0.5, 0.5}});
                        mag_results.idx = arma::urowvec(4, arma::fill::ones);
                        mag_results.S = arma::rowvec(4, arma::fill::ones);
                        break;
                    case 2:
                        atom_results.r = arma::mat({{0, 0.5, 0,   0.5, 0.1, 0.6, 0.9, 0.4, 0.9, 0.4, 0.1, 0.6},
                                                    {0, 0.5, 0,   0.5, 0.1, 0.6, 0.9, 0.4, 0.1, 0.6, 0.9, 0.4},
                                                    {0, 0,   0.5, 0.5, 0,   0,   0.5, 0.5, 0.5, 0.5, 0,   0}});
                        atom_results.idx = arma::urowvec({1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2});
                        atom_results.mag = arma::urowvec(12, arma::fill::ones);
                        mag_results.r = arma::mat({{0, 0.5, 0,   0.5, 0.1, 0.6, 0.9, 0.4, 0.9, 0.4, 0.1, 0.6},
                                                   {0, 0.5, 0,   0.5, 0.1, 0.6, 0.9, 0.4, 0.1, 0.6, 0.9, 0.4},
                                                   {0, 0,   0.5, 0.5, 0,   0,   0.5, 0.5, 0.5, 0.5, 0,   0}});
                        mag_results.idx = arma::urowvec({1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2});
                        mag_results.S = arma::rowvec({1, 1, 1, 1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5});
                        break;
                    case 3:
                        atom_results.r = arma::mat({{0, 0.5, 0,   0.5, 0.1, 0.6, 0.9, 0.4, 0.9, 0.4, 0.1, 0.6},
                                                    {0, 0.5, 0,   0.5, 0.1, 0.6, 0.9, 0.4, 0.1, 0.6, 0.9, 0.4},
                                                    {0, 0,   0.5, 0.5, 0,   0,   0.5, 0.5, 0.5, 0.5, 0,   0}});
                        atom_results.idx = arma::urowvec({1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2});
                        atom_results.mag = arma::urowvec({0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1});
                        mag_results.r = arma::mat({{0.1, 0.6, 0.9, 0.4, 0.9, 0.4, 0.1, 0.6},
                                                   {0.1, 0.6, 0.9, 0.4, 0.1, 0.6, 0.9, 0.4},
                                                   {0,   0,   0.5, 0.5, 0.5, 0.5, 0,   0}});
                        mag_results.idx = 2*arma::urowvec(8, arma::fill::ones);
                        mag_results.S = 0.5*arma::rowvec(8, arma::fill::ones);
                        break;
                }
                break;

            case 84:
                switch (noAtoms) {
                    case 1:
                        atom_results.r = arma::mat({{0, 0.5, 0.5, 0},
                                                    {0, 0.5, 0,   0.5},
                                                    {0, 0,   0,   0}});
                        atom_results.idx = arma::urowvec(4, arma::fill::ones);
                        atom_results.mag = arma::urowvec(4, arma::fill::ones);
                        mag_results.r = arma::mat({{0, 0.5, 0.5, 0},
                                                   {0, 0.5, 0,   0.5},
                                                   {0, 0,   0,   0}});
                        mag_results.idx = arma::urowvec(4, arma::fill::ones);
                        mag_results.S = arma::rowvec(4, arma::fill::ones);
                        break;
                    case 2:
                        atom_results.r = arma::mat({{0, 0.5, 0.5, 0,   0.1, 0.4, 0.4, 0.1, 0.9, 0.6, 0.6, 0.9},
                                                    {0, 0.5, 0,   0.5, 0.1, 0.4, 0.1, 0.4, 0.9, 0.6, 0.9, 0.6},
                                                    {0, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0}});
                        atom_results.idx = arma::urowvec({1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2});
                        atom_results.mag = arma::urowvec(12, arma::fill::ones);
                        mag_results.r = arma::mat({{0, 0.5, 0.5, 0,   0.1, 0.4, 0.4, 0.1, 0.9, 0.6, 0.6, 0.9},
                                                   {0, 0.5, 0,   0.5, 0.1, 0.4, 0.1, 0.4, 0.9, 0.6, 0.9, 0.6},
                                                   {0, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0}});
                        mag_results.idx = arma::urowvec({1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2});
                        mag_results.S = arma::rowvec({ 1, 1, 1, 1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 });
                        break;
                    case 3:
                        atom_results.r = arma::mat({{0, 0.5, 0.5, 0,   0.1, 0.4, 0.4, 0.1, 0.9, 0.6, 0.6, 0.9},
                                                    {0, 0.5, 0,   0.5, 0.1, 0.4, 0.1, 0.4, 0.9, 0.6, 0.9, 0.6},
                                                    {0, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0}});
                        atom_results.idx = arma::urowvec({1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2});
                        atom_results.mag = arma::urowvec({ 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1 });
                        mag_results.r = arma::mat({{0.1, 0.4, 0.4, 0.1, 0.9, 0.6, 0.6, 0.9},
                                                   {0.1, 0.4, 0.1, 0.4, 0.9, 0.6, 0.9, 0.6},
                                                   {0,   0,   0,   0,   0,   0,   0,   0}});
                        mag_results.idx = arma::urowvec({2, 2, 2, 2, 2, 2, 2, 2});
                        mag_results.S = 0.5*arma::rowvec(8, arma::fill::ones);
                        break;
                }
                break;
        }
    }
};
//// A new one of these is create for each test
class SwCoreAtom : public testing::TestWithParam<SwCoreAtomResults> {
public:
    virtual void SetUp() {}
    virtual void TearDown() {}
};

INSTANTIATE_TEST_CASE_P(OneMagAtom,
                        SwCoreAtom,
                        ::testing::Values(SwCoreAtomResults(0, 1),
                                          SwCoreAtomResults(9, 1),
                                          SwCoreAtomResults(19, 1),
                                          SwCoreAtomResults(84, 1)
                        )
);

INSTANTIATE_TEST_CASE_P(TwoMagAtom,
                        SwCoreAtom,
                        ::testing::Values(SwCoreAtomResults(0, 2),
                                          SwCoreAtomResults(9, 2),
                                          SwCoreAtomResults(19, 2),
                                          SwCoreAtomResults(84, 2)
                        )
);
INSTANTIATE_TEST_CASE_P(OneMagOneNot,
                        SwCoreAtom,
                        ::testing::Values(SwCoreAtomResults(0, 3),
                                          SwCoreAtomResults(9, 3),
                                          SwCoreAtomResults(19, 3),
                                          SwCoreAtomResults(84, 3)
                        )
);
//
TEST_P(SwCoreAtom, ATOM){
    SwCoreAtomResults theseResults = GetParam();

    unit_cell thisUnitCell = unit_cell();
    lattice thisLattice = lattice();

    thisUnitCell.r[0] = 0;
    thisUnitCell.r[1] = 0;
    thisUnitCell.r[2] = 0;
    thisUnitCell.nAtom = 1;
    thisUnitCell.S[0] = 1;

    switch (theseResults.noAtoms) {
        case 1:
            break;
        case 2:
            thisUnitCell.r[3] = 0.1;
            thisUnitCell.r[4] = 0.1;
            thisUnitCell.r[5] = 0;
            thisUnitCell.nAtom = 2;
            thisUnitCell.S[1] = 0.5;
            break;
        case 3:
            thisUnitCell.r[3] = 0.1;
            thisUnitCell.r[4] = 0.1;
            thisUnitCell.r[5] = 0;
            thisUnitCell.nAtom = 2;
            thisUnitCell.S[0] = 0;
            thisUnitCell.S[1] = 0.5;
            break;
    }

    thisLattice.nSymOp = theseResults.symIDX;
    spinw s = spinw(thisLattice, thisUnitCell, twin(), matrix(), single_ion(), coupling(), mag_str(), unit());


    std::tuple<arma::mat, arma::urowvec, arma::urowvec> resultTP = s.atom();
    atom_struct expectedTP = theseResults.atom_results;

    arma::mat A  = expectedTP.r;
    arma::mat AA = get<0>(resultTP);

    std::vector<double> ExpectedA(A.memptr(), A.memptr() + A.n_elem);
    std::vector<double> ResultA(AA.memptr(), AA.memptr() + AA.n_elem);

    arma::urowvec B  = expectedTP.idx;
    arma::urowvec BB = get<1>(resultTP);
    std::vector<arma::uword> ExpectedB(B.memptr(), B.memptr() + B.n_elem);
    std::vector<arma::uword> ResultB(BB.memptr(), BB.memptr() + BB.n_elem);

    arma::urowvec C  = expectedTP.mag;
    arma::urowvec CC = get<2>(resultTP);
    std::vector<arma::uword> ExpectedC(C.memptr(), C.memptr() + C.n_elem);
    std::vector<arma::uword> ResultC(CC.memptr(), CC.memptr() + CC.n_elem);

    EXPECT_EQ(ExpectedA, ResultA);
    EXPECT_EQ(ExpectedB, ResultB);
    EXPECT_EQ(ExpectedC, ResultC);
}

TEST_P(SwCoreAtom, MATOM) {
    SwCoreAtomResults theseResults = GetParam();

    unit_cell thisUnitCell = unit_cell();
    lattice thisLattice = lattice();

    thisUnitCell.r[0] = 0;
    thisUnitCell.r[1] = 0;
    thisUnitCell.r[2] = 0;
    thisUnitCell.nAtom = 1;
    thisUnitCell.S[0] = 1;

    switch (theseResults.noAtoms) {
        case 1:
            break;
        case 2:
            thisUnitCell.r[3] = 0.1;
            thisUnitCell.r[4] = 0.1;
            thisUnitCell.r[5] = 0;
            thisUnitCell.nAtom = 2;
            thisUnitCell.S[1] = 0.5;
            break;
        case 3:
            thisUnitCell.r[3] = 0.1;
            thisUnitCell.r[4] = 0.1;
            thisUnitCell.r[5] = 0;
            thisUnitCell.nAtom = 2;
            thisUnitCell.S[0] = 0;
            thisUnitCell.S[1] = 0.5;
            break;
    }

    thisLattice.nSymOp = theseResults.symIDX;
    spinw s = spinw(thisLattice, thisUnitCell, twin(), matrix(), single_ion(), coupling(), mag_str(), unit());

    matom_struct resultMA = s.matom();
    matom_struct expectedMA = theseResults.mag_results;

    arma::mat D  = expectedMA.r;
    arma::mat DD = resultMA.r;
    std::vector<double> ExpectedD(D.memptr(), D.memptr() + D.n_elem);
    std::vector<double> ResultD(DD.memptr(), DD.memptr() + DD.n_elem);

    arma::urowvec E = expectedMA.idx;
    arma::urowvec EE = resultMA.idx;
    std::vector<arma::uword > ExpectedE(E.memptr(), E.memptr() + E.n_elem);
    std::vector<arma::uword> ResultE(EE.memptr(), EE.memptr() + EE.n_elem);

    arma::rowvec F = expectedMA.S;
    arma::rowvec FF = resultMA.S;
    std::vector<double> ExpectedF(F.memptr(), F.memptr() + F.n_elem);
    std::vector<double> ResultF(FF.memptr(), FF.memptr() + FF.n_elem);

    EXPECT_EQ(ExpectedD, ResultD);
    EXPECT_EQ(ExpectedE, ResultE);
    EXPECT_EQ(ExpectedF, ResultF);
}