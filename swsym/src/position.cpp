#include <../include/position.h>

void position(TYPE symOp, TYPE r0, TYPE fid, TYPE tol, TYPE& r, TYPE& aIdx, _Opinfo& opInfo)
{
  TYPE idx, idx2, idxF, ii, isMoved, nAtom, nGenAtom, opMove, rTemp ;
  rowvec _aux_rowvec_1 ;
  if (m2cpp::isempty(symOp))
  {
    symOp = arma::join_rows(arma::eye<mat>(3, 3), arma::zeros<mat>(3, 1)) ;
  }
  if (size(r0, 1)!=3)
  {
    std::cerr << "position:WrongInput"<< "The positions have to be in column vector format!" << std::endl ;
  }
  nAtom = size(r0, 2) ;
  r = arma::zeros<mat>(3, 0) ;
  aIdx.reset() ;
  isMoved = cell(1, nAtom) ;
  opMove = arma::zeros<cube>(3, 3, 0) ;
  for (ii=1; ii<=nAtom; ii++)
  {
    double __aux_rowvec_1 [] = {1, 3, 2} ;
    _aux_rowvec_1 = rowvec(__aux_rowvec_1, 3, false) ;
    rTemp = permute(mod(mmat(symOp(m2cpp::span<uvec>(0, symOp.n_rows-1), m2cpp::fspan(1, 1, 3), m2cpp::span<uvec>(0, symOp.n_slices-1)), r0.cols(ii))+symOp(m2cpp::span<uvec>(0, symOp.n_rows-1), 4, m2cpp::span<uvec>(0, symOp.n_slices-1)), 1), _aux_rowvec_1) ;
    isMoved{} = arma::sum(arma::square(bsxfun(@minus, sw_cmod(rTemp, tol), sw_cmod(r0.cols(ii), tol))), 1)>pow(tol, 2) ;
    if (numel(rTemp)>3)
    {
      [rTemp, idxF] = sw_uniquetol(sw_cmod(rTemp, tol), tol) ;
    }
    else
    {
      idxF = 1 ;
    }
    r = {arma::join_rows(r, rTemp)} ;
    aIdx = {arma::join_rows(aIdx, arma::ones<rowvec>(size(rTemp, 2))*ii)} ;
    opMove = cat(3, opMove, symOp(m2cpp::span<uvec>(0, symOp.n_rows-1), m2cpp::fspan(1, 1, 3), idxF)) ;
  }
  opInfo.ismoved = isMoved ;
  opInfo.opmove = opMove ;
  nGenAtom = numel(aIdx) ;
  if (fid!=0)
  {
    idx = 1 ;
    for (ii=1; ii<=nAtom; ii++)
    {
      std::printf(fid, "\nAtomic coordinates generated for: (%5.3f %5.3f %5.3f)\n", r0.cols(ii)) ;
      idx2 = 1 ;
      while (idx<=nGenAtom && aIdx(idx)==ii)
      {
        std::printf(fid, "R%02i = (%5.3f %5.3f %5.3f)\n", idx2, r.cols(idx)) ;
        idx2 = idx2+1 ;
        idx = idx+1 ;
      }
    }
  }
}

TYPE sw_cmod(TYPE r, TYPE tol)
{
  
  r = mod(r, 1) ;
  r(r>1-tol) = r(r>1-tol)-1 ;
  return r ;
}