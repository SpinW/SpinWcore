//
// Created by ward_s on 11/04/18.
//

#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "../../include/templateFuncs.h"

//using ::testing::Return;
//using ::testing::_;
using namespace std;

TEST(TemplateTest, MMAT1)
{
    arma::mat A = {1, 2, 3, 4, 5};
    arma::mat B = A;
    A = repmat(A.t(),3,1);

    int dim[] = {1, 2};

    arma::mat C = {
             {1.0000,    2.0000,    3.0000,    4.0000,    5.0000},
             {2.0000,    4.0000,    6.0000,    8.0000,   10.0000},
             {3.0000,    6.0000,    9.0000,   12.0000,   15.0000},
             {4.0000,    8.0000,   12.0000,   16.0000,   20.0000},
             {5.0000,   10.0000,   15.0000,   20.0000,   25.0000},
             {1.0000,    2.0000,    3.0000,    4.0000,    5.0000},
             {2.0000,    4.0000,    6.0000,    8.0000,   10.0000},
             {3.0000,    6.0000,    9.0000,   12.0000,   15.0000},
             {4.0000,    8.0000,   12.0000,   16.0000,   20.0000},
             {5.0000,   10.0000,   15.0000,   20.0000,   25.0000},
             {1.0000,    2.0000,    3.0000,    4.0000,    5.0000},
             {2.0000,    4.0000,    6.0000,    8.0000,   10.0000},
             {3.0000,    6.0000,    9.0000,   12.0000,   15.0000},
             {4.0000,    8.0000,   12.0000,   16.0000,   20.0000},
             {5.0000,   10.0000,   15.0000,   20.0000,   25.0000}};


    arma::mat CC = mmat(A,B,dim);
    auto *Expected = new double[C.n_elem];
    memcpy(Expected, C.memptr(), sizeof(double) * C.n_elem);
    auto *Result = new double[CC.n_elem];
    memcpy(Result, CC.memptr(), sizeof(double) * CC.n_elem);
    EXPECT_EQ(*Expected,*Result);
}

TEST(TemplateTest, MMAT2)
{
    arma::mat A = {{1, 2, 3, 4, 5}};
    A = repmat(A,4,1);
    arma::mat B = {{1, 2, 3, 4}};
    int dim[] = {2, 1};

    arma::mat C = {{10,    20,    30,    40,    50}};
    arma::mat CC = mmat(A,B,dim);

    auto *Expected = new double[C.n_elem];
    memcpy(Expected, C.memptr(), sizeof(double) * C.n_elem);
    auto *Result = new double[CC.n_elem];
    memcpy(Result, CC.memptr(), sizeof(double) * CC.n_elem);
    EXPECT_EQ(*Expected,*Result);
}

TEST(TemplateTest, MMAT3)
{
    arma::mat A = {{1, 2, 3, 4, 5}};
    A = repmat(A,4,1);
    arma::mat B = {{1, 2, 3, 4}};
    int dim[] = {1, 2};

    arma::mat C = {10,    20,    30,    40,    50};
    C = C.t();
    A = A.t();
    B = B.t();

    arma::mat CC = mmat(A,B,dim);

    auto *Expected = new double[C.n_elem];
    memcpy(Expected, C.memptr(), sizeof(double) * C.n_elem);
    auto *Result = new double[CC.n_elem];
    memcpy(Result, CC.memptr(), sizeof(double) * CC.n_elem);
    EXPECT_EQ(*Expected,*Result);
}

TEST(TemplateTest, MMAT4)
{
    arma::mat A = {{1, 2, 3, 4}};
    arma::mat B = A.t();
    int dim[] = {1, 2};

    arma::mat C = {30};

    arma::mat CC = mmat(A,B,dim);

    auto *Expected = new double[C.n_elem];
    memcpy(Expected, C.memptr(), sizeof(double) * C.n_elem);
    auto *Result = new double[CC.n_elem];
    memcpy(Result, CC.memptr(), sizeof(double) * CC.n_elem);
    EXPECT_EQ(*Expected,*Result);
}

TEST(TemplateTest, MMAT5)
{
    arma::mat A = {{1, 2, 3, 4}};
    arma::mat B = A.t();
    int dim[] = {2, 1};

    arma::mat C = {{1,     2,     3,     4},
                   {2,     4,     6,     8},
                   {3,     6,     9,    12},
                   {4,     8,    12,    16}};

    arma::mat CC = mmat(A,B,dim);

    auto *Expected = new double[C.n_elem];
    memcpy(Expected, C.memptr(), sizeof(double) * C.n_elem);
    auto *Result = new double[CC.n_elem];
    memcpy(Result, CC.memptr(), sizeof(double) * CC.n_elem);
    EXPECT_EQ(*Expected,*Result);
}