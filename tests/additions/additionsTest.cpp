//
// Created by ward_s on 11/04/18.
//
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "../../include/sw_additions.h"

using namespace std;
#define PI 3.14159265

TEST(SW_Structs, latticeTest){
    lattice* lattice1 = new lattice();
    double angle[] = {90, 90, 90};
    double lat_const[] = {3.3, 4.4, 6};

    for (int i = 0; i < 3; i++) {
        lattice1->angle[i] = angle[i];
        EXPECT_DOUBLE_EQ(angle[i], lattice1->angle[i]);
        lattice1->lat_const[i] = lat_const[i];
        EXPECT_DOUBLE_EQ(lat_const[i], lattice1->lat_const[i]);
    }
}

class SW_AdditionsClass : public testing::TestWithParam<bool> {
public:
    virtual void SetUp() {}

    virtual void TearDown() {}
};

INSTANTIATE_TEST_CASE_P(SW_Additions,
                        SW_AdditionsClass,
                        ::testing::Values(true, false)
);


TEST_P(SW_AdditionsClass, basisVector){

    lattice lattice1 = lattice();
    double angle[] = {90, 90, 120};
    double lat_const[] = {3.3, 4.4, 6};

    for (int i = 0; i < 3; i++) {
        lattice1.angle[i] = 2*PI*angle[i]/360;
        lattice1.lat_const[i] = lat_const[i];
    }

    bool option = GetParam();
    arma::mat Result_C  = arma_basisvector(lattice1,option);

    arma::mat Expected_C(3,3);
    if (option){
        Expected_C = {{1.0000, -0.5000, 0.0000},
                      {0,       0.8660, 0.0000},
                      {0,       0,      1}};
    } else {
        Expected_C = {{3.3000, -2.2000, 0.0000},
                      {0,       3.8105, 0.0000},
                      {0,       0,      6}};
    }

    std::vector<double> Result(Result_C.memptr(),
                               Result_C.memptr() + Result_C.n_elem);

    std::vector<double> Expected(Expected_C.memptr(),
                                 Expected_C.memptr() + Expected_C.n_elem);


    ASSERT_EQ(Expected.size(),Result.size());
    for (int i = 0; i < Result.size(); i++){
        EXPECT_NEAR(Result[i], Expected[i],1E-4); // Because maths is hard :-/
    }
}
