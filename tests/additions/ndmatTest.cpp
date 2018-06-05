//
// Created by ward_s on 04/06/18.
//
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "../../include/ndmat.h"

using namespace std;
#define PI 3.14159265

//
// Created by ward_s on 29/05/18.
//
//
//#include <armadillo>
//#include "ndmat.h"
//
//int main(int argc, char** argv)
//{
//    arma::vec inVec = arma::linspace(1, 32, 32);
//    int dims[] = {2, 2, 2, 1, 2, 2};
//    int D = 6;
//
//    ndmat in = ndmat(inVec.memptr(), &(dims[0]), D);
//
//    std::cout << "Checking input is correct (MATLAB formulism)...." << std::endl;
//
//    arma::mat vis1 = inVec.t();
//    vis1.insert_rows(vis1.n_rows, arma::conv_to<arma::mat>::from((in.getAll()).t()));
//    vis1.insert_rows(vis1.n_rows, arma::conv_to<arma::rowvec>::from((in.getIndex() + 1).t()));
//    vis1.insert_rows(vis1.n_rows, arma::conv_to<arma::mat>::from(in.getSubs() + 1));
//
//    std::cout << vis1 << std::endl;
//
//
//    std::cout << "Doing a perumation" << std::endl;
//    int perm[] = {1, 0, 2, 3, 4, 5};
//    in.permute(perm);
//    arma::mat vis2 = inVec.t();
//    vis2.insert_rows(vis2.n_rows, arma::conv_to<arma::mat>::from((in.getAll()).t()));
//    vis2.insert_rows(vis2.n_rows, arma::conv_to<arma::rowvec>::from((in.getIndex() + 1).t()));
//    vis2.insert_rows(vis2.n_rows, arma::conv_to<arma::mat>::from(in.getSubs() + 1));
//
//    std::cout << vis2 << std::endl;
//
//    std::cout << "Doing a repmat" << std::endl;
//    ndmat in2 = ndmat(inVec.memptr(), &(dims[0]), D);
//    int expand[] = {1, 1, 2, 1, 1, 1};
//    in2.repmat(expand);
//    arma::mat vis3 = arma::conv_to<arma::mat>::from((in2.getAll()).t());
//    vis3.insert_rows(vis3.n_rows, arma::conv_to<arma::rowvec>::from((in2.getIndex() + 1).t()));
//    vis3.insert_rows(vis3.n_rows, arma::conv_to<arma::mat>::from(in2.getSubs() + 1));
//
//    std::cout << vis3 << std::endl;
//
//    std::cout << "Removing some of a dimension" << std::endl;
//    int p = 2;
//    int rmRange[] = {2, 3}; // As it was 0, 1
//    in2.remove1D(p, rmRange);
//    arma::mat vis4 = arma::conv_to<arma::mat>::from((in2.getAll()).t());
//    vis4.insert_rows(vis4.n_rows, arma::conv_to<arma::rowvec>::from((in2.getIndex() + 1).t()));
//    vis4.insert_rows(vis4.n_rows, arma::conv_to<arma::mat>::from(in2.getSubs() + 1));
//
//    std::cout << vis4 << std::endl;
//
//    std::cout << "Reshaping the tensor" << std::endl;
//    int nNewDims = 2;
//    int newDims[] = {2, 16};
//    in2.reshape(nNewDims, newDims);
//    arma::mat vis5 = arma::conv_to<arma::mat>::from((in2.getAll()).t());
//    vis5.insert_rows(vis5.n_rows, arma::conv_to<arma::rowvec>::from((in2.getIndex() + 1).t()));
//    vis5.insert_rows(vis5.n_rows, arma::conv_to<arma::mat>::from(in2.getSubs() + 1));
//
//    std::cout << vis5 << std::endl;
//
//    std::cout << "Multiplying" << std::endl;
//    nNewDims = 3;
//    int newDims2[] = {2, 4, 4};
//    in2.reshape(nNewDims, newDims2); // Make it a better size to play around with.
//
//    arma::vec mVec = arma::linspace(1, 4, 4);
//    int mDims[] = {1, 4, 1};
//    int mD = 3;
//    ndmat m1 = ndmat(mVec.memptr(), mDims, mD);
//    in2.multiply(m1);
//    arma::mat vis6 = arma::conv_to<arma::mat>::from((in2.getAll()).t());
//    vis6.insert_rows(vis6.n_rows, arma::conv_to<arma::rowvec>::from((in2.getIndex() + 1).t()));
//    vis6.insert_rows(vis6.n_rows, arma::conv_to<arma::mat>::from(in2.getSubs() + 1));
//
//    std::cout << vis6 << std::endl;
//
//
//    std::cout << "Summation, dim 1 (MATLAB 2)" << std::endl;
//    in2.sum(1);
//
//    arma::mat vis7 = arma::conv_to<arma::mat>::from((in2.getAll()).t());
//    vis7.insert_rows(vis7.n_rows, arma::conv_to<arma::rowvec>::from((in2.getIndex() + 1).t()));
//    vis7.insert_rows(vis7.n_rows, arma::conv_to<arma::mat>::from(in2.getSubs() + 1));
//
//    std::cout << vis7 << std::endl;
//
//    arma::vec inVec1 = arma::linspace(0, 0.31, 32);
//
//    ndcx_mat inI = ndcx_mat(inVec.memptr(), inVec1.memptr(), &(dims[0]), D);
//
//    std::cout << "Doing a perumation IMAG" << std::endl;
//    inI.permute(perm);
//    std::cout << "Imag_Matrix: \n" << (inI.getAll()).t() << std::endl;
//    arma::mat vis8 = arma::conv_to<arma::mat>::from(inVec.t());
//    vis8.insert_rows(vis8.n_rows, arma::conv_to<arma::mat>::from(inVec1.t()));
//    vis8.insert_rows(vis8.n_rows, arma::conv_to<arma::rowvec>::from((inI.getIndex() + 1).t()));
//    vis8.insert_rows(vis8.n_rows, arma::conv_to<arma::mat>::from(inI.getSubs() + 1));
//
//    std::cout << "Index Data: \n" << vis8 << std::endl;
//    return 0;
//}
