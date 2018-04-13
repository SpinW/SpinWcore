//
// Created by ward_s on 11/04/18.
//
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "../../include/templateFuncs.h"

using namespace std;

TEST(TemplateTest, MMAT1)
{
    arma::mat A = {1, 2, 3, 4, 5};
    arma::mat B = A;
    A = repmat(A.t(), 3, 1);

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
            {5.0000,   10.0000,   15.0000,   20.0000,   25.0000}
    };

    arma::mat CC = mmat(A,B,dim);

    std::vector<double> Expected(C.memptr(), C.memptr() + C.n_elem);
    std::vector<double> Result(CC.memptr(), CC.memptr() + CC.n_elem);
    EXPECT_EQ(Expected, Result);
}

TEST(TemplateTest, MMAT2)
{
    arma::mat A = {{1, 2, 3, 4, 5}};
    A = repmat(A, 4, 1);
    arma::mat B = {{1, 2, 3, 4}};
    int dim[] = {2, 1};

    arma::mat C = {{10, 20, 30, 40, 50}};
    arma::mat CC = mmat(A, B, dim);

    std::vector<double> Expected(C.memptr(), C.memptr() + C.n_elem);
    std::vector<double> Result(CC.memptr(), CC.memptr() + CC.n_elem);

    EXPECT_EQ(Expected, Result);
}

TEST(TemplateTest, MMAT3)
{
    arma::mat A = {{1, 2, 3, 4, 5}};
    A = repmat(A, 4, 1);
    arma::mat B = {{1, 2, 3, 4}};
    int dim[] = {1, 2};

    arma::mat C = {10, 20, 30, 40, 50};

    C = C.t();
    A = A.t();
    B = B.t();

    arma::mat CC = mmat(A, B, dim);

    std::vector<double> Expected(C.memptr(), C.memptr() + C.n_elem);
    std::vector<double> Result(CC.memptr(), CC.memptr() + CC.n_elem);
    EXPECT_EQ(Expected, Result);
}

TEST(TemplateTest, MMAT4)
{
    arma::mat A = {{1, 2, 3, 4}};
    arma::mat B = A.t();
    int dim[] = {1, 2};

    arma::mat C = {30};

    arma::mat CC = mmat(A, B, dim);

    std::vector<double> Expected(C.memptr(), C.memptr() + C.n_elem);
    std::vector<double> Result(CC.memptr(), CC.memptr() + CC.n_elem);
    EXPECT_EQ(Expected, Result);
}

TEST(TemplateTest, MMAT5)
{
    arma::mat A = {{1, 2, 3, 4}};
    arma::mat B = A.t();
    int dim[] = {2, 1};

    arma::mat C = {
            {1, 2, 3,  4 },
            {2, 4, 6,  8 },
            {3, 6, 9,  12},
            {4, 8, 12, 16}
    };

    arma::mat CC = mmat(A, B, dim);

    std::vector<double> Expected(C.memptr(), C.memptr() + C.n_elem);
    std::vector<double> Result(CC.memptr(), CC.memptr() + CC.n_elem);
    EXPECT_EQ(Expected, Result);
}

TEST(TemplateTest, MMAT_DIM1)
{
    arma::mat A = {{1, 2, 3, 4}};
    arma::mat B = A.t();
    int dim[] = {1, 1};

    arma::mat C = {30};

    EXPECT_THROW(
            {
                try {
                    arma::mat CC = mmat(A, B, dim);
                } catch (const std::invalid_argument &e) {
                    EXPECT_STREQ("Dimension element is duplicated", e.what());
                    throw;
                }
            }, std::invalid_argument);
}

TEST(TemplateTest, MMAT_DIM2)
{
    arma::mat A = {{1, 2, 3, 4}};
    arma::mat B = A.t();
    int dim[] = {1, 3};

    EXPECT_THROW(
            {
                try {
                    arma::mat CC = mmat(A, B, dim);
                } catch (const std::invalid_argument &e) {
                    EXPECT_STREQ("Dimension element has to be 1 or 2", e.what());
                    throw;
                }
            }, std::invalid_argument);

    dim[0] = 3; dim[1] = 1;
    EXPECT_THROW(
            {
                try {
                    arma::mat CC = mmat(A, B, dim);
                } catch (const std::invalid_argument &e) {
                    EXPECT_STREQ("Dimension element has to be 1 or 2", e.what());
                    throw;
                }
            }, std::invalid_argument);
    dim[0] = -1; dim[1] = 1;
    EXPECT_THROW(
            {
                try {
                    arma::mat CC = mmat(A, B, dim);
                } catch (const std::invalid_argument &e) {
                    EXPECT_STREQ("Dimension element has to be 1 or 2", e.what());
                    throw;
                }
            }, std::invalid_argument);

    dim[0] = 1; dim[1] = -1;
    EXPECT_THROW(
            {
                try {
                    arma::mat CC = mmat(A, B, dim);
                } catch (const std::invalid_argument &e) {
                    EXPECT_STREQ("Dimension element has to be 1 or 2", e.what());
                    throw;
                }
            }, std::invalid_argument);

    dim[0] = -1; dim[1] = -2;
    EXPECT_THROW(
            {
                try {
                    arma::mat CC = mmat(A, B, dim);
                } catch (const std::invalid_argument &e) {
                    EXPECT_STREQ("Dimension element has to be 1 or 2", e.what());
                    throw;
                }
            }, std::invalid_argument);
}

TEST(TemplateTest, UNIQUETOL1)
{
    using namespace arma;
    arma::mat A = randu<mat>(3, 10);
    double delta = 0.02;
    A(span(), 4) = A(span(), 0) + 0.25*delta;
    A(span(), 6) = A(span(), 3) + 0.1*delta;

    arma::mat C = A;
    C.shed_col(4);
    C.shed_col(5);
    arma::mat CC = uniquetol(A, delta);

    std::vector<double> Expected(C.memptr(), C.memptr() + C.n_elem);
    std::vector<double> Result(CC.memptr(), CC.memptr() + CC.n_elem);
    EXPECT_EQ(Expected, Result);
}

arma::cube createTestCube(arma::uword i, arma::uword j,arma::uword k) {
    arma::uword l = 0;
    arma::cube A(i, j, k, arma::fill::ones);
    for (int ii = 0; ii < i; ii++) {
        for (int jj = 0; jj < j; jj++) {
            for (int kk = 0; kk < k; kk++) {
                A(ii, jj, kk) = l;
                l++;
            }
        }
    }
    return A;
}


TEST(TemplateTest, PERMUTE_123) {
    using namespace arma;

    uword i, j, k;
    i = 2; j = 3; k = 1;

    cube A = createTestCube(i, j, k);

    int opt[3] = {1, 2, 3};
    auto B = permute(A, opt);

    std::vector<arma::uword> C_sz  = {B.n_rows, B.n_cols, B.n_slices};
    std::vector<arma::uword> CC_sz = {i, j, k};

    EXPECT_EQ(C_sz, CC_sz);

    arma::cube C(CC_sz[0], CC_sz[1], CC_sz[2], arma::fill::zeros);
    C.slice(0) = {
            {0, 1, 2},
            {3, 4, 5}
    };

    std::vector<double> Result(B.memptr(),   B.memptr() + B.n_elem);
    std::vector<double> Expected(C.memptr(), C.memptr() + C.n_elem);

    EXPECT_EQ(Expected, Result);
}

TEST(TemplateTest, PERMUTE_132) {
    using namespace arma;

    uword i, j, k;
    i = 2; j = 3; k = 1;

    cube A = createTestCube(i, j, k);

    int opt[3] = {1, 3, 2};
    auto B = permute(A, opt);

    std::vector<arma::uword> C_sz  = {B.n_rows, B.n_cols, B.n_slices};
    std::vector<arma::uword> CC_sz = {i, k, j};

    EXPECT_EQ(C_sz,CC_sz);

    arma::cube C(CC_sz[0], CC_sz[1], CC_sz[2], arma::fill::zeros);
    C.slice(0) = arma::colvec({0, 3});
    C.slice(1) = arma::colvec({1, 4});
    C.slice(2) = arma::colvec({2, 5});

    std::vector<double> Result(B.memptr(),   B.memptr() + B.n_elem);
    std::vector<double> Expected(C.memptr(), C.memptr() + C.n_elem);

    EXPECT_EQ(Expected, Result);
}

TEST(TemplateTest, PERMUTE_213) {
    using namespace arma;

    uword i, j, k;
    i = 2; j = 3; k = 1;

    cube A = createTestCube(i, j, k);

    int opt[3] = {2, 1, 3};
    auto B = permute(A, opt);

    std::vector<arma::uword> C_sz  = {B.n_rows, B.n_cols, B.n_slices};
    std::vector<arma::uword> CC_sz = {j, i, k};

    EXPECT_EQ(C_sz,CC_sz);

    arma::cube C(CC_sz[0], CC_sz[1], CC_sz[2], arma::fill::zeros);
    C.slice(0) = {
            {0, 3},
            {1, 4},
            {2, 5}
    };

    std::vector<double> Result(B.memptr(),   B.memptr() + B.n_elem);
    std::vector<double> Expected(C.memptr(), C.memptr() + C.n_elem);

    EXPECT_EQ(Expected, Result);
}

TEST(TemplateTest, PERMUTE_231) {
    using namespace arma;

    uword i, j, k;
    i = 2; j = 3; k = 1;

    cube A = createTestCube(i, j, k);

    int opt[3] = {2, 3, 1};
    auto B = permute(A, opt);

    std::vector<arma::uword> C_sz  = {B.n_rows, B.n_cols, B.n_slices};
    std::vector<arma::uword> CC_sz = {j, k, i};

    EXPECT_EQ(C_sz,CC_sz);

    arma::cube C(CC_sz[0], CC_sz[1], CC_sz[2], arma::fill::zeros);
    C.slice(0) = arma::colvec({0, 1, 2});
    C.slice(1) = arma::colvec({3, 4, 5});

    std::vector<double> Result(B.memptr(),   B.memptr() + B.n_elem);
    std::vector<double> Expected(C.memptr(), C.memptr() + C.n_elem);

    EXPECT_EQ(Expected, Result);
}

TEST(TemplateTest, PERMUTE_321) {
    using namespace arma;

    uword i, j, k;
    i = 2; j = 3; k = 1;

    cube A = createTestCube(i, j, k);

    int opt[3] = {3, 2, 1};
    auto B = permute(A, opt);

    std::vector<arma::uword> C_sz  = {B.n_rows, B.n_cols, B.n_slices};
    std::vector<arma::uword> CC_sz = {k, j, i};

    EXPECT_EQ(C_sz,CC_sz);

    arma::cube C(CC_sz[0], CC_sz[1], CC_sz[2], arma::fill::zeros);
    C.slice(0) = {0, 1, 2};
    C.slice(1) = {3, 4, 5};

    std::vector<double> Result(B.memptr(),   B.memptr() + B.n_elem);
    std::vector<double> Expected(C.memptr(), C.memptr() + C.n_elem);

    EXPECT_EQ(Expected, Result);
}

TEST(TemplateTest, PERMUTE_312) {
    using namespace arma;

    uword i, j, k;
    i = 2; j = 3; k = 1;

    cube A = createTestCube(i, j, k);

    int opt[3] = {3, 1, 2};
    auto B = permute(A, opt);

    std::vector<arma::uword> C_sz  = {B.n_rows, B.n_cols, B.n_slices};
    std::vector<arma::uword> CC_sz = {k, i, j};

    EXPECT_EQ(C_sz,CC_sz);

    arma::cube C(CC_sz[0], CC_sz[1], CC_sz[2], arma::fill::zeros);
    C.slice(0) = {{0, 3}};
    C.slice(1) = {{1, 4}};
    C.slice(2) = {{2, 5}};

    std::vector<double> Result(B.memptr(),   B.memptr() + B.n_elem);
    std::vector<double> Expected(C.memptr(), C.memptr() + C.n_elem);

    EXPECT_EQ(Expected, Result);
}