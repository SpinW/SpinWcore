//
// Created by ward_s on 11/04/18.
//
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "../../swsym/include/swsym.h"

using namespace std;
const char actualDir[] = "/MATLAB/mtools/SpinW_Dev/spinw/dat_files/";

class symResults{
public:
    arma::cube Expected;
    int Index;
    arma::cube Result(){
        return thisSym.getSym(Index);
    }
    explicit symResults(int ind){
        switch (ind){
            case 0:{
                arma::cube temp(3,4,1,arma::fill::zeros);
                temp.slice(0) = {{1, 0, 0, 0},
                                     {0, 1, 0, 0},
                                     {0, 0, 1, 0}};
                Expected = temp;}
                break;
            case 9: {
                arma::cube temp(3, 4, 2, arma::fill::zeros);
                temp.slice(0) = {{-1, 0, 0,  0},
                                 {0,  1, 0,  0},
                                 {0,  0, -1, 0}};

                temp.slice(1) = {{-1, 0,  0,  0},
                                     {0,  -1, 0,  0},
                                     {0,  0,  -1, 0}};
                Expected = temp;
            }
                break;
            case 19: {
                arma::cube temp(3, 4, 3, arma::fill::zeros);

                temp.slice(0) = {{1, 0, 0, 0.5},
                                     {0, 1, 0, 0.5},
                                     {0, 0, 1, 0.0}};

                temp.slice(1) = {{-1, 0,  0, 0.0},
                                     {0,  -1, 0, 0.0},
                                     {0,  0,  1, 0.5}};

                temp.slice(2) = {{-1, 0, 0,  0.0},
                                     {0,  1, 0,  0.0},
                                     {0,  0, -1, 0.5}};
                Expected = temp;
            }
                break;
            case 39: {
                arma::cube temp(3, 4, 3, arma::fill::zeros);

                temp.slice(0) = {{1, 0, 0, 0.0},
                                     {0, 1, 0, 0.5},
                                     {0, 0, 1, 0.5}};

                temp.slice(1) = {{-1, 0,  0, 0.0},
                                     {0,  -1, 0, 0.0},
                                     {0,  0,  1, 0.0}};

                temp.slice(2) = {{    1, 0,  0,  0.5},
                                     {0, -1, 0,  0.0},
                                     {0, 0,  1, 0.0}};
                Expected = temp;
            }
                break;
            case 79: {
                arma::cube temp(3, 4, 3, arma::fill::zeros);

                temp.slice(0) = {{1, 0, 0, 0.5},
                                     {0, 1, 0, 0.5},
                                     {0, 0, 1, 0.5}};

                temp.slice(1) = {{-1, 0,  0, 0.5},
                                     {0,  -1, 0, 0.5},
                                     {0,  0,  1, 0.5}};

                temp.slice(2) = {{0, -1, 0, 0.0},
                                     {1, 0,  0, 0.5},
                                     {0, 0,  1, 0.25}};
                Expected = temp;
            }
                break;
            case 159: {
                arma::cube temp(3, 4, 3, arma::fill::zeros);

                temp.slice(0) = {{1, 0, 0, 1.0 / 3},
                                     {0, 1, 0, 2.0 / 3},
                                     {0, 0, 1, 2.0 / 3}};

                temp.slice(1) = {{0, -1, 0, 0.0},
                                     {1, -1, 0, 0.0},
                                     {0, 0,  1, 0.0}};

                temp.slice(2) = {{0,  -1, 0, 0.0},
                                     {-1, 0,  0, 0.0},
                                     {0,  0,  1, 0.0}};
                Expected = temp;
            }
                break;
        }
        Index = ind;
    }
protected:
    char tempAdr[256];
    swsym thisSym = swsym(strcpy(tempAdr,string(actualDir).c_str()));
};

// A new one of these is create for each test
class SwSymTestClassLoc : public testing::TestWithParam<const char*>
{
public:
    virtual void SetUp(){}
    virtual void TearDown(){}
};

class SwSymTestClassSym : public testing::TestWithParam<symResults>
{
public:
    virtual void SetUp(){}
    virtual void TearDown(){}
};

INSTANTIATE_TEST_CASE_P(DatDirs,
                        SwSymTestClassLoc,
                        ::testing::Values(actualDir,
                                          "/MATLAB/mtools/SpinW_Dev/spinw",
                                          "/test/fail/"
                        ));

INSTANTIATE_TEST_CASE_P(SymPos,
                        SwSymTestClassSym,
                        ::testing::Values(symResults(0),
                                          symResults(9),
                                          symResults(19),
                                          symResults(39),
                                          symResults(79),
                                          symResults(159)
                        ));

TEST_P(SwSymTestClassLoc, acceptsDir)
{
    const char *dir = GetParam();
    string temp = string(dir);
    char blah[256];
    try {
        swsym* thisSym = new swsym(strcpy(blah,temp.c_str()));
    } catch (const std::exception &e){
        cout << temp << endl;
        EXPECT_THROW({
            EXPECT_STREQ("Unable to open dat file",e.what());
            throw;
        },std::exception);
    }
}

TEST_P(SwSymTestClassSym, testSymResult) {

    symResults thisSymResults = GetParam();

    arma::cube Result_C   = thisSymResults.Result();
    arma::cube Expected_C = thisSymResults.Expected;

    std::vector<double> Result(Result_C.memptr(),
                               Result_C.memptr() + Result_C.n_elem);

    std::vector<double> Expected(Expected_C.memptr(),
                                 Expected_C.memptr() + Expected_C.n_elem);
//    for (size_t ind = 0; ind < Result.size(); ind++){
//        EXPECT_DOUBLE_EQ(Result[ind],Expected[ind]);
//    }
    EXPECT_EQ(Result,Expected);
}