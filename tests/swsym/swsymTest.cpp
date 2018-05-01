//
// Created by ward_s on 11/04/18.
//
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "../../swsym/include/swsym.h"

#pragma clang diagnostic push
#pragma ide diagnostic ignored "TemplateArgumentsIssues"
using namespace std;

class symResults {
public:
    int Index = -1;
//    swsym thisSym = swsym(strcpy(tempAdr, string(actualDir).c_str()));
    swsym thisSym = swsym();

    explicit symResults(int ind) {

        if (ind >= 0){
            Index = ind;
        } else {
            throw std::range_error(std::string("Index must be greater than or equal to 0"));
        }
    }

    arma::cube Result() {
        if (Index >= 0){
            return thisSym.getSym(Index);
        } else {
            throw std::range_error(std::string("Index must be greater than or equal to 0"));
        }
    }

    arma::cube Expected(){
        arma::cube temp;
        switch (Index) {
            case 0:
                temp = arma::cube(3, 4, 1);
                temp.slice(0) = {{1, 0, 0, 0},
                                 {0, 1, 0, 0},
                                 {0, 0, 1, 0}};
                break;
            case 9:
                temp = arma::cube(3, 4, 2);
                temp.slice(0) = {{-1, 0, 0,  0},
                                 {0,  1, 0,  0},
                                 {0,  0, -1, 0}};

                temp.slice(1) = {{-1, 0,  0,  0},
                                 {0,  -1, 0,  0},
                                 {0,  0,  -1, 0}};
                break;
            case 19:
                temp = arma::cube(3, 4, 3);
                temp.slice(0) = {{1, 0, 0, 0.5},
                                 {0, 1, 0, 0.5},
                                 {0, 0, 1, 0.0}};

                temp.slice(1) = {{-1, 0,  0, 0.0},
                                 {0,  -1, 0, 0.0},
                                 {0,  0,  1, 0.5}};

                temp.slice(2) = {{-1, 0, 0,  0.0},
                                 {0,  1, 0,  0.0},
                                 {0,  0, -1, 0.5}};
                break;
            case 39:
                temp = arma::cube(3, 4, 3);
                temp.slice(0) = {{1, 0, 0, 0.0},
                                 {0, 1, 0, 0.5},
                                 {0, 0, 1, 0.5}};

                temp.slice(1) = {{-1, 0,  0, 0.0},
                                 {0,  -1, 0, 0.0},
                                 {0,  0,  1, 0.0}};

                temp.slice(2) = {{1, 0,  0, 0.5},
                                 {0, -1, 0, 0.0},
                                 {0, 0,  1, 0.0}};
                break;
            case 69:
                temp = arma::cube(3, 4, 5);
                temp.slice(0) = {{1, 0, 0, 0.5},
                                 {0, 1, 0, 0.5},
                                 {0, 0, 1, 0.0}};

                temp.slice(1) = {{1, 0, 0, 0.5},
                                 {0, 1, 0, 0.0},
                                 {0, 0, 1, 0.5}};

                temp.slice(2) = {{-1,  0, 0, 0.75},
                                 { 0, -1, 0, 0.75},
                                 { 0,  0, 1, 0.0}};

                temp.slice(3) = {{-1, 0,  0, 0.75},
                                 { 0, 1,  0, 0.0},
                                 { 0, 0, -1, 0.75}};

                temp.slice(4) = {{-1,  0,  0, 0},
                                 { 0, -1,  0, 0},
                                 { 0,  0, -1, 0}};
                break;
            case 79:
                temp = arma::cube(3, 4, 3);
                temp.slice(0) = {{1, 0, 0, 0.5},
                                 {0, 1, 0, 0.5},
                                 {0, 0, 1, 0.5}};

                temp.slice(1) = {{-1, 0,  0, 0.5},
                                 {0,  -1, 0, 0.5},
                                 {0,  0,  1, 0.5}};

                temp.slice(2) = {{0, -1, 0, 0.0},
                                 {1, 0,  0, 0.5},
                                 {0, 0,  1, 0.25}};
                break;
            case 97:
                temp = arma::cube(3, 4, 4);
                temp.slice(0) = {{1, 0, 0, 0.5},
                                 {0, 1, 0, 0.5},
                                 {0, 0, 1, 0.5}};

                temp.slice(1) = {{-1,  0, 0, 0.5},
                                 { 0, -1, 0, 0.5},
                                 { 0,  0, 1, 0.5}};

                temp.slice(2) = {{0, -1, 0, 0.0},
                                 {1,  0, 0, 0.5},
                                 {0,  0, 1, 0.25}};

                temp.slice(3) = {{-1, 0,  0, 0.5},
                                 { 0, 1,  0, 0.0},
                                 { 0, 0, -1, 0.75}};
                break;
            case 140:
                temp = arma::cube(3, 4, 5);
                temp.slice(0) = {{1, 0, 0, 0.5},
                                 {0, 1, 0, 0.5},
                                 {0, 0, 1, 0.5}};

                temp.slice(1) = {{-1,  0, 0, 0.5},
                                 { 0, -1, 0, 0.0},
                                 { 0,  0, 1, 0.5}};

                temp.slice(2) = {{0, -1, 0, 0.25},
                                 {1,  0, 0, 0.75},
                                 {0,  0, 1, 0.25}};

                temp.slice(3) = {{-1, 0,  0, 0.5},
                                 { 0, 1,  0, 0.0},
                                 { 0, 0, -1, 0.5}};

                temp.slice(4) = {{-1,  0,  0, 0},
                                 { 0, -1,  0, 0},
                                 { 0,  0, -1, 0}};
                break;
            case 159:
                temp = arma::cube(3, 4, 3);
                temp.slice(0) = {{1, 0, 0, 1.0 / 3},
                                 {0, 1, 0, 2.0 / 3},
                                 {0, 0, 1, 2.0 / 3}};

                temp.slice(1) = {{0, -1, 0, 0.0},
                                 {1, -1, 0, 0.0},
                                 {0, 0,  1, 0.0}};

                temp.slice(2) = {{0,  -1, 0, 0.0},
                                 {-1, 0,  0, 0.0},
                                 {0,  0,  1, 0.0}};
                break;
        }
        return temp;
    }


protected:
    char tempAdr[256];
};


class symPosResults {
public:
    int Index = -1;
    swsym thisSym = swsym();
    double toll = 1E-4;

    explicit symPosResults(int ind) {

        if (ind >= 0){
            Index = ind;
        } else {
            throw std::range_error(std::string("Index must be greater than or equal to 0"));
        }
    }

    std::tuple<arma::mat, arma::urowvec> Result(arma::mat r0) {
        if (Index >= 0){
            return this->thisSym.position(r0,Index,toll);
        } else {
            throw std::range_error(std::string("Index must be greater than or equal to 0"));
        }
    }

    std::tuple<arma::mat, arma::Row<int>> Expected(){
        arma::mat temp;
        arma::Row<int> idx;
        switch (Index) {
            case 14:
                temp = arma::mat({{0.5000, 0, 0, 0.6000, 0.9000, 0.9000},
                                  {0.5000, 0, 0, 0.6000, 0.1000, 0.9000},
                                  {0, 0.5000, 0, 0, 0.5000, 0}});
                idx = arma::Row<int>({0, 0, 0, 1, 1, 1});
                break;
            case 2:
                temp = arma::mat({{0, 0.9000},
                                  {0, 0.1000},
                                  {0, 0}});
                idx = arma::Row<int>({0, 1});
                break;
            case 89:
                temp = arma::mat({{0, 0.5, 0.5, 0.9, 0.4, 0.4},
                                  {0, 0.5, 0.5, 0.9, 0.6, 0.6},
                                  {0, 0, 0, 0, 0, 0}});
                idx = arma::Row<int>({0, 0, 0, 1, 1, 1});
                break;
            case 4:
                temp = arma::mat({{0.5, 0, 0.6, 0.9},
                                  {0.5, 0, 0.6, 0.1},
                                  {0, 0, 0, 0}});
                idx = arma::Row<int>({0, 0, 1, 1});
                break;
        }
        return make_tuple(temp,idx);
    }


protected:
    char tempAdr[256];
};

// A new one of these is create for each test
class SwSymTestClassLoc : public testing::TestWithParam<const char *> {
public:
    virtual void SetUp() {}

    virtual void TearDown() {}
};

class SwSymTestClassSym : public testing::TestWithParam<symResults> {
public:
    virtual void SetUp() {}

    virtual void TearDown() {}
};

class SwSymTestClass : public testing::TestWithParam<const char *> {
public:
    virtual void SetUp() {}

    virtual void TearDown() {}
};
class SwSymTestClassPos : public testing::TestWithParam<int> {
public:
    virtual void SetUp() {}

    virtual void TearDown() {}
};

//INSTANTIATE_TEST_CASE_P(SwSym,
//                        SwSymTestClassLoc,
//                        ::testing::Values(actualDir,
//                                          "/MATLAB/mtools/SpinW_Dev/spinw",
//                                          "/test/fail/"
//                        )
//);

INSTANTIATE_TEST_CASE_P(SwSym,
                        SwSymTestClassSym,
                        ::testing::Values(symResults(0),
                                          symResults(9),
                                          symResults(19),
                                          symResults(39),
                                          symResults(79),
                                          symResults(159)
                        )
);

INSTANTIATE_TEST_CASE_P(SwSym,
                        SwSymTestClassPos,
                        ::testing::Values(14, 2, 89, 4)
);
INSTANTIATE_TEST_CASE_P(SwSym,
                        SwSymTestClass,
                        ::testing::Values(
                                "069x+1/2,y+1/2,z; x+1/2,y,z+1/2; -x+3/4,-y+3/4,z; -x+3/4,y,-z+3/4; -x,-y,-z",
                                "097x+1/2,y+1/2,z+1/2; -x+1/2,-y+1/2,z+1/2; -y,x+1/2,z+1/4; -x+1/2,y,-z+3/4",
                                "140x+1/2,y+1/2,z+1/2; -x+1/2,-y,z+1/2; -y+1/4,x+3/4,z+1/4; -x+1/2,y,-z+1/2; -x,-y,-z"
                        )
);

//TEST_P(SwSymTestClassLoc, testDatLoading) {
//    const char *dir = GetParam();
//    string temp = string(dir);
//    char blah[256];
//    try {
//        swsym thisSym = swsym(strcpy(blah, temp.c_str()));
//    } catch (const std::exception &e) {
//        if (temp == string(actualDir)){
//            EXPECT_THROW({throw;},std::exception);
//        } else {
//            EXPECT_THROW({
//                             EXPECT_STREQ("Unable to open dat file", e.what());
//                             throw;
//                         }, std::exception);
//        }
//    }
//}

TEST_P(SwSymTestClassSym, testDatInterpreting) {

    symResults thisSymResults = GetParam();

    arma::cube Result_C = thisSymResults.Result();
    arma::cube Expected_C = thisSymResults.Expected();

    std::vector<double> Result(Result_C.memptr(),
                               Result_C.memptr() + Result_C.n_elem);

    std::vector<double> Expected(Expected_C.memptr(),
                                 Expected_C.memptr() + Expected_C.n_elem);

    ASSERT_EQ(Expected.size(),Result.size());
    EXPECT_EQ(Result, Expected);
}

TEST_P(SwSymTestClass, testStrRead) {
    const char *symChar = GetParam();
    string temp = string(symChar);
    string thisID = temp.substr(0,3);
    string thisSymStr = temp.substr(3,temp.length());

    int ID;
    sscanf(thisID.c_str(), "%d", &ID);

    auto thisResult = symResults(ID);
    arma::cube Expected_C(3,4,300,arma::fill::zeros);
    thisResult.thisSym.interpretSymString(Expected_C,thisSymStr);

    arma::cube Result_C = thisResult.Result();

    std::vector<double> Result(Result_C.memptr(),
                               Result_C.memptr() + Result_C.n_elem);

    std::vector<double> Expected(Expected_C.memptr(),
                                 Expected_C.memptr() + Expected_C.n_elem);

    ASSERT_EQ(Expected.size(),Result.size());
    EXPECT_EQ(Result, Expected);
}
TEST_P(SwSymTestClass, testStrAdd) {
    const char *symChar = GetParam();
    string temp = string(symChar);
    string thisID = temp.substr(0,3);
    string thisSymStr = temp.substr(3,temp.length());

    int ID;
    sscanf(thisID.c_str(), "%d", &ID);

    auto thisResult = symResults(ID);

    thisResult.thisSym.addSymString(thisSymStr);
    arma::cube Result_C = thisResult.thisSym.getSym(thisResult.thisSym.totalSyms-1);
    arma::cube Expected_C = thisResult.Expected();

    std::vector<double> Result(Result_C.memptr(),
                               Result_C.memptr() + Result_C.n_elem);

    std::vector<double> Expected(Expected_C.memptr(),
                                 Expected_C.memptr() + Expected_C.n_elem);

    ASSERT_EQ(Expected.size(),Result.size());
    EXPECT_EQ(Result, Expected);
}


TEST_P(SwSymTestClassPos,Position){
    int ID = GetParam();
    arma::mat r0(3,2,arma::fill::zeros);
    r0(arma::span(),1) = arma::colvec({0.1, 0.1, 0.0});
    double toll = 1E-4;
    auto thisResult = symPosResults(ID);
    thisResult.toll = toll;
    std::tuple<arma::mat, arma::Row<int>> Expected_C = thisResult.Expected();
    std::tuple<arma::mat, arma::urowvec> Result_C = thisResult.Result(r0);

    arma::mat R1 = get<0>(Result_C);
    arma::urowvec R2 = get<1>(Result_C);

    arma::mat E1 = get<0>(Expected_C);
    arma::Row<int> E2 = get<1>(Expected_C);

    ASSERT_EQ(R1.size(),E1.size());
    ASSERT_EQ(R2.size(),E2.size());

    std::vector<double> resultPos(R1.memptr(), R1.memptr() + R1.n_elem);
    std::vector<double> expectedPos(E1.memptr(), E1.memptr() + E1.n_elem);

    std::vector<int> resultIdx(R2.memptr(), R2.memptr() + R2.n_elem);
    std::vector<int> expectedIdx(E2.memptr(), E2.memptr() + E2.n_elem);

    EXPECT_EQ(resultPos, expectedPos);
    EXPECT_EQ(resultIdx, expectedIdx);
}

TEST(SwSymTestBond, Bond){
    arma::mat R  = arma::mat({{0.5, 0, 0.5}, {0, 0.5, 0.5},  {0, 0, 0}});
    arma::mat BV = arma::mat({{6, -3, 0},    {0, 5.1962, 0}, {0, 0, 5}});

    arma::mat E1 = { {0,  0, 0, 0, 1, -1},
                     {1, -1, 0, 0, 0,  0},
                     {0,  0, 0, 0, 0,  0},
                     {2,  0, 1, 2, 0,  1},
                     {0,  1, 2, 0, 1,  2}};

    arma::urowvec E2 = {1, 1, 1, 1, 1, 1};

// bond is composed of [dl_x, dl_y, dl_z, start_idx, end_idx, number];
    arma::colvec bond = arma::colvec({0, 1, 0, 2, 0, 1});

    swsym thisSym = swsym();
    int index = thisSym.searchSym(string("P -3"));
    double toll = 1E-5;

    std::tuple<arma::mat, arma::urowvec> Results =  thisSym.bond(R, BV, bond, index, toll);

    arma::mat R1 = std::get<0>(Results);
    arma::urowvec R2 = std::get<1>(Results);

    std::vector<double> resultBond(R1.memptr(), R1.memptr() + R1.n_elem);
    std::vector<double> expectedBond(E1.memptr(), E1.memptr() + E1.n_elem);

    std::vector<arma::uword> resultIdx(R2.memptr(), R2.memptr() + R2.n_elem);
    std::vector<arma::uword> expectedIdx(E2.memptr(), E2.memptr() + E2.n_elem);

    ASSERT_EQ(R1.size(),E1.size());
    ASSERT_EQ(R2.size(),E2.size());

    EXPECT_EQ(resultBond, expectedBond);
    EXPECT_EQ(resultIdx, expectedIdx);

}

#pragma clang diagnostic pop
#pragma clang diagnostic pop