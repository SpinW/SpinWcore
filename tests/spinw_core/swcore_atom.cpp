//
// Created by ward_s on 03/05/18.
//

#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "../../src/templateFuncs.tpp"
#include "../../src/spinw.cpp"

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

    SwCoreAtomResults(int inSym) {
        symIDX = inSym;
        switch (inSym) {
            case 0:
                atom_results.r = arma::mat(3, 1, arma::fill::zeros);
                atom_results.idx = arma::urowvec(1, arma::fill::zeros);
                atom_results.mag = arma::urowvec(1, arma::fill::ones);
                mag_results.r = arma::mat(3, 1, arma::fill::zeros);
                mag_results.idx = arma::urowvec(1, arma::fill::zeros);
                mag_results.S = arma::rowvec(1, arma::fill::ones);
                break;
            case 9:
                atom_results.r = arma::mat(3, 1, arma::fill::zeros);
                atom_results.idx = arma::urowvec(1, arma::fill::zeros);
                atom_results.mag = arma::urowvec(1, arma::fill::ones);
                mag_results.r = arma::mat(3, 1, arma::fill::zeros);
                mag_results.idx = arma::urowvec(1, arma::fill::zeros);
                mag_results.S = arma::rowvec(1, arma::fill::ones);
                break;

            case 19:
                atom_results.r = arma::mat({{0.5, 0},{0.5, 0},{0, 0.5}});
                atom_results.idx = arma::urowvec(2, arma::fill::zeros);
                atom_results.mag = arma::urowvec(2, arma::fill::ones);
                mag_results.r = arma::mat({{0.5, 0},
                                           {0.5, 0},
                                           {0, 0.5}});
                mag_results.idx = arma::urowvec(2, arma::fill::zeros);
                mag_results.S = arma::rowvec(2, arma::fill::ones);
                break;
            case 84:
                atom_results.r = arma::mat({{0.5, 0.5, 0},{0.5, 0, 0},{0, 0, 0}});
                atom_results.idx = arma::urowvec(3, arma::fill::zeros);
                atom_results.mag = arma::urowvec(3, arma::fill::ones);
                mag_results.r = arma::mat({{0.5, 0.5, 0},
                                           {0.5, 0, 0},
                                           {0, 0, 0}});
                mag_results.idx = arma::urowvec(3, arma::fill::zeros);
                mag_results.S = arma::rowvec(3, arma::fill::ones);
                break;
        }
    }
};
//// A new one of these is create for each test
class SwCoreAtom : public testing::TestWithParam<SwCoreAtomResults> {
public:
    virtual void SetUp() {}
    virtual void TearDown() {}
//    static void SetUpTestCase() {
//        thisUnitCell = new unit_cell();
//        thisLattice = new lattice();
//
//        thisUnitCell->r[0][0] = 0;
//        thisUnitCell->r[1][0] = 0;
//        thisUnitCell->r[2][0] = 0;
//        thisUnitCell->nAtom = 1;
//        thisUnitCell->S[0] = 1;
//        thisLattice->nSymOp = -1;
//    }
//    static void TearDownTestCase() {
//        delete thisLattice;
//        thisLattice = NULL;
//        delete thisUnitCell;
//        thisUnitCell = NULL;
//    }
//
//    static lattice* thisLattice;
//    static unit_cell* thisUnitCell;
};
//
//lattice* SwCoreAtom::thisLattice = NULL;
//unit_cell* SwCoreAtom::thisUnitCell = NULL;
//
INSTANTIATE_TEST_CASE_P(SwCore,
                        SwCoreAtom,
                        ::testing::Values(SwCoreAtomResults(0),
                                          SwCoreAtomResults(9),
                                          SwCoreAtomResults(19),
                                          SwCoreAtomResults(84)
                        )
);
//
TEST_P(SwCoreAtom, ATOM){
    SwCoreAtomResults theseResults = GetParam();

    unit_cell thisUnitCell = unit_cell();
    lattice thisLattice = lattice();

    thisUnitCell.r[0][0] = 0;
    thisUnitCell.r[1][0] = 0;
    thisUnitCell.r[2][0] = 0;
    thisUnitCell.nAtom = 1;
    thisUnitCell.S[0] = 1;
    thisLattice.nSymOp = theseResults.symIDX;
    spinw s = spinw(thisLattice, thisUnitCell, twin(), mag_str(), unit());


    std::tuple<arma::mat, arma::urowvec, arma::urowvec> resultTP = s.atom();
    atom_struct expectedTP = theseResults.atom_results;

    arma::mat A  = expectedTP.r;
    arma::mat AA = get<0>(resultTP);

    std::cout << AA << std::endl;

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

TEST_P(SwCoreAtom, MATOM){
    SwCoreAtomResults theseResults = GetParam();

    unit_cell thisUnitCell = unit_cell();
    lattice thisLattice = lattice();

    thisUnitCell.r[0][0] = 0;
    thisUnitCell.r[1][0] = 0;
    thisUnitCell.r[2][0] = 0;
    thisUnitCell.nAtom = 1;
    thisUnitCell.S[0] = 1;
    thisLattice.nSymOp = theseResults.symIDX;
    spinw s = spinw(thisLattice, thisUnitCell, twin(), mag_str(), unit());

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
    std::vector<arma::uword > ExpectedF(F.memptr(), F.memptr() + F.n_elem);
    std::vector<arma::uword> ResultF(FF.memptr(), FF.memptr() + FF.n_elem);

    EXPECT_EQ(ExpectedD, ResultD);
    EXPECT_EQ(ExpectedE, ResultE);
    EXPECT_EQ(ExpectedF, ResultF);
}

//TEST(SwCoreAtom, ATOM) {
//
//    lattice thisLattice = lattice();
//    unit_cell thisUnitCell = unit_cell();
//    thisUnitCell.r[0][0] = 0;
//    thisUnitCell.r[1][0] = 0;
//    thisUnitCell.r[2][0] = 0;
//    thisUnitCell.nAtom = 1;
//    thisUnitCell.S[0] = 1;
//
//    thisLattice.nSymOp = 0;
//
//
//    spinw s = spinw(thisLattice, thisUnitCell, twin(), mag_str(), unit());
//    std::tuple<arma::mat, arma::urowvec, arma::urowvec> resultTP = s.atom();
//
//    arma::mat AA = get<0>(resultTP);
//    std::vector<double> Expected(3, 0.0);
//    std::vector<double> Result(AA.memptr(), AA.memptr() + AA.n_elem);
//
//    arma::urowvec BB = get<1>(resultTP);
//    arma::urowvec CC = get<2>(resultTP);
//
//    EXPECT_EQ(Expected, Result);
//    EXPECT_EQ(BB(0), 0);
//    EXPECT_EQ(CC(0), 1);
//}
//
//TEST(SwCoreAtom, MATOM) {
//    lattice thisLattice = lattice();
//    unit_cell thisUnitCell = unit_cell();
//    thisUnitCell.r[0][0] = 0;
//    thisUnitCell.r[1][0] = 0;
//    thisUnitCell.r[2][0] = 0;
//    thisUnitCell.nAtom = 1;
//    thisUnitCell.S[0] = 1;
//
//    thisLattice.nSymOp = 0;
//
//
//    spinw s = spinw(thisLattice, thisUnitCell, twin(), mag_str(), unit());
//    matom_struct resultTP = s.matom();
//
//    arma::mat AA = resultTP.r;
//    std::vector<double> Expected(3, 0.0);
//    std::vector<double> Result(AA.memptr(), AA.memptr() + AA.n_elem);
//
//    arma::urowvec BB = resultTP.idx;
//    arma::rowvec CC = resultTP.S;
//
//    EXPECT_EQ(Expected, Result);
//    EXPECT_EQ(BB(0), 0);
//    EXPECT_EQ(CC(0), 1);
//}