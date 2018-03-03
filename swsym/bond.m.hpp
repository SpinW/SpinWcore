// Automatically translated using m2cpp 2.0 on 2018-03-01 12:10:35

#ifndef BOND_M_HPP
#define BOND_M_HPP

#include "mconvert.h"
#include <iostream>
#include <armadillo>
using namespace arma ;

void bond(TYPE r, TYPE bv, TYPE bond, TYPE symOp, TYPE tol, TYPE& genCp, TYPE& ugenCp) ;
void isnewUC(TYPE A, TYPE B, TYPE tol, TYPE& isnew, TYPE& symIdx) ;
TYPE sw_cmod(TYPE r, TYPE tol) ;
TYPE cfloor(TYPE r0, TYPE tol) ;
TYPE uniqueb(TYPE bond) ;

void bond(TYPE r, TYPE bv, TYPE bond, TYPE symOp, TYPE tol, TYPE& genCp, TYPE& ugenCp)
{
  TYPE atom1, atom2, dist, dl, dlnew, iNew, r1, r1new, r2, r2new, rightDist, tolDist ;
  rowvec _aux_rowvec_1, _aux_rowvec_2, _aux_rowvec_3 ;
  tolDist = 1e-5 ;
  r1 = r.cols(bond(4)) ;
  r2 = r.cols(bond(5)) ;
  bond(m2cpp::fspan(1, 1, 3), dl) ;
  double __aux_rowvec_1 [] = {1, 3, 2} ;
  _aux_rowvec_1 = rowvec(__aux_rowvec_1, 3, false) ;
  r1new = permute(mmat(symOp(m2cpp::span<uvec>(0, symOp.n_rows-1), m2cpp::fspan(1, 1, 3), m2cpp::span<uvec>(0, symOp.n_slices-1)), r1)+symOp(m2cpp::span<uvec>(0, symOp.n_rows-1), 4, m2cpp::span<uvec>(0, symOp.n_slices-1)), _aux_rowvec_1) ;
  double __aux_rowvec_2 [] = {1, 3, 2} ;
  _aux_rowvec_2 = rowvec(__aux_rowvec_2, 3, false) ;
  r2new = permute(mmat(symOp(m2cpp::span<uvec>(0, symOp.n_rows-1), m2cpp::fspan(1, 1, 3), m2cpp::span<uvec>(0, symOp.n_slices-1)), r2)+symOp(m2cpp::span<uvec>(0, symOp.n_rows-1), 4, m2cpp::span<uvec>(0, symOp.n_slices-1)), _aux_rowvec_2) ;
  double __aux_rowvec_3 [] = {1, 3, 2} ;
  _aux_rowvec_3 = rowvec(__aux_rowvec_3, 3, false) ;
  dlnew = permute(mmat(symOp(m2cpp::span<uvec>(0, symOp.n_rows-1), m2cpp::fspan(1, 1, 3), m2cpp::span<uvec>(0, symOp.n_slices-1)), dl), _aux_rowvec_3)-cfloor(r1new, tol)+cfloor(r2new, tol) ;
  r1new = mod(r1new, 1) ;
  r2new = mod(r2new, 1) ;
  isnewUC(r, r1new, tolDist, iNew, atom1) ;
  if (any(iNew))
  {
    std::cerr << "bond:SymProblem"<< "The generated positions for atom1 are wrong!" << std::endl ;
  }
  isnewUC(r, r2new, tolDist, iNew, atom2) ;
  if (any(iNew))
  {
    std::cerr << "bond:SymProblem"<< "The generated positions for atom2 are wrong!" << std::endl ;
  }
  dist = sqrt(arma::sum(arma::square((bv*(r.cols(atom2)-r.cols(atom1)+dlnew))), 1)) ;
  rightDist = abs(dist-dist(1))<tol ;
  if (!all(rightDist))
  {
    warning("bond:SymProblem", "Symmetry generated couplings are dropped!") ;
  }
  genCp = {dlnew, atom1, atom2} ;
  genCp = genCp.cols(rightDist) ;
  if (nargout>1)
  {
    ugenCp = uniqueb(genCp) ;
  }
}

void isnewUC(TYPE A, TYPE B, TYPE tol, TYPE& isnew, TYPE& symIdx)
{
  TYPE idx, nA, nB, notequal ;
  rowvec _aux_rowvec_1, _aux_rowvec_2 ;
  nA = size(A, 2) ;
  nB = size(B, 2) ;
  double __aux_rowvec_1 [] = {2, 3, 1} ;
  _aux_rowvec_1 = rowvec(__aux_rowvec_1, 3, false) ;
  double __aux_rowvec_2 [] = {3, 2, 1} ;
  _aux_rowvec_2 = rowvec(__aux_rowvec_2, 3, false) ;
  notequal = arma::sum(arma::square(sw_cmod(abs(repmat(permute(A, _aux_rowvec_1), {arma::join_rows(arma::join_rows(m2cpp::srow<sword>(1), nB), m2cpp::srow<sword>(1))})-repmat(permute(B, _aux_rowvec_2), {arma::join_rows(arma::join_rows(nA, m2cpp::srow<sword>(1)), m2cpp::srow<sword>(1))})), tol)), 3)>tol ;
  isnew = all(notequal, 1) ;
  idx = m2cpp::fspan(1, 1, nB) ;
  symIdx = max(bsxfun(@times, !notequal.cols(idx(!isnew)), arma::trans((m2cpp::fspan(1, 1, size(notequal, 1))))), {}, 1) ;
}

TYPE sw_cmod(TYPE r, TYPE tol)
{
  
  r = mod(r, 1) ;
  r(r>1-tol) = r(r>1-tol)-1 ;
  return r ;
}

TYPE cfloor(TYPE r0, TYPE tol)
{
  TYPE idx, r ;
  r = arma::floor(r0) ;
  idx = abs(r0-r)>1-tol ;
  r(idx) = r(idx)+1 ;
  return r ;
}

TYPE uniqueb(TYPE bond)
{
  TYPE c1, c2, nC, nc1, uniqueB ;
  rowvec _aux_rowvec_1, _aux_rowvec_2, _aux_rowvec_3, _aux_rowvec_4 ;
  nC = size(bond(), 2) ;
  double __aux_rowvec_1 [] = {2, 3, 1} ;
  _aux_rowvec_1 = rowvec(__aux_rowvec_1, 3, false) ;
  c1 = permute(bond(), _aux_rowvec_1) ;
  double __aux_rowvec_2 [] = {3, 2, 1} ;
  _aux_rowvec_2 = rowvec(__aux_rowvec_2, 3, false) ;
  c2 = permute(bond(), _aux_rowvec_2) ;
  double __aux_rowvec_3 [] = {5, 4} ;
  _aux_rowvec_3 = rowvec(__aux_rowvec_3, 2, false) ;
  double __aux_rowvec_4 [] = {2, 3, 1} ;
  _aux_rowvec_4 = rowvec(__aux_rowvec_4, 3, false) ;
  nc1 = permute({-bond(m2cpp::fspan(1, 1, 3), m2cpp::span<uvec>(0, bond.n_cols-1)), bond(_aux_rowvec_3, m2cpp::span<uvec>(0, bond.n_cols-1))}, _aux_rowvec_4) ;
  uniqueB = all(trimatu(any(bsxfun(@ne, c1, c2), 3) && any(bsxfun(@ne, nc1, c2), 3)) || trimatl(arma::ones<umat>(nC)), 1) ;
  return uniqueB ;
}
#endif