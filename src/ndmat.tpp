//
// Created by ward_s on 28/05/18.
//
#include <armadillo>

template <typename T, arma::uword aSize = 1>
class ndMat {

public:
    arma::umat subs;
    arma::uvec inds;
    arma::uvec D;
    arma::Mat<typename T::elem_type> dataStore;

public:
    // Dummy to avoid hassle. Init with ndMat<arma::___>();
    explicit ndMat(){
        int dims[] = {1};
        ndMat(arma::Mat<typename T::elem_type>(aSize, 1),&(dims[0]) , 1);
    };

    ndMat(arma::Mat<typename T::elem_type> inMat, int *dims, int nDims){
        D = arma::conv_to<arma::uvec>::from(arma::Col<int>(&(dims[0]), nDims));
        dataStore = (inMat.n_cols > 1) ? arma::reshape(inMat, inMat.n_elem, 1) : inMat;
        inds = arma::linspace<arma::uvec>(0, dataStore.n_elem - 1, dataStore.n_elem);
        subs = ind2sub(inds);
    }

    template <typename TT>
    void multiply(TT byMat){
        arma::Col<int>rMat = arma::conv_to<arma::Col<int>>::from(D - byMat.D + 1);
        byMat.repmat(rMat.memptr());
        dataStore %= byMat.dataStore;
    };

    template <typename TT>
    void divide(TT byMat){
        arma::Col<int>rMat = arma::conv_to<arma::Col<int>>::from(D - byMat.D + 1);
        byMat.repmat(rMat.memptr());
        dataStore /= byMat.dataStore;
    };


    void repmat(int* dims){
        arma::uvec thisD = arma::conv_to<arma::uvec>::from(arma::Col<int>(&(dims[0]), D.n_rows));
        arma::uword TF = arma::prod(thisD);
        dataStore = arma::repmat(dataStore, TF, 1);
        for (arma::uword i = 0; i < D.n_elem; i++){
            if (thisD(i) > 1){
                arma::umat sub_base = subs;
                for (arma::uword j = 1; j < thisD(i); j++){
                    arma::uvec pf(D.n_elem, arma::fill::zeros);
                    pf(i) = j*D(i);
                    sub_base.each_col() += pf;
                    subs.insert_cols(subs.n_cols, sub_base);
                }
            }
        }
        D %= thisD;
        inds = sub2ind(subs);
        if (!inds.is_sorted()) {
            arma::uvec sortI = arma::sort_index(inds);
            inds = inds.elem(sortI);
            subs = subs.cols(sortI);
            dataStore = dataStore.rows(sortI);
        }
    }

    void permute(int *dims) {
        arma::uvec thisIdx = arma::conv_to<arma::uvec>::from(arma::Col<int>(&(dims[0]), D.n_rows));
        subs = subs.rows(thisIdx);
        inds = sub2ind(subs);
    }

    void remove(int numDim, int* dimensions, int* range){
        for (int i = 0; i < numDim; i++){
            // We set a remove dimension to 1, so dimensions input is not changed.
            remove1D(dimensions[i], &(range[i*2]));
        }
    }

    void remove1D(int dimension, int* range){
        arma::urowvec mySubs = subs.row(dimension);
        arma::urowvec R = (mySubs < range[0]) || (mySubs > range[1]);
        dataStore = dataStore.rows(arma::find(R));
        inds = arma::linspace<arma::uvec>(0, dataStore.n_elem -1, dataStore.n_elem);
        D(dimension) -= (range[1] - range[0] + 1);
        subs = ind2sub(inds);
    }

    void reshape(int nDims, int* dims){
        arma::Col<int>temp(&(dims[0]), nDims);
        arma::uvec tempI = temp < 0;
        if (arma::any(tempI)){
            arma::uvec fInd = arma::find(tempI);
            if (fInd.n_elem == 1) {
                temp.row(fInd(0)) = arma::prod(D) / arma::prod(D.elem(arma::find(-1 * (tempI - 1))));
            } else {
                std::cout << "You made a bo bo" << std::endl;
            }
        }
        D = arma::conv_to<arma::uvec>::from(temp);
        subs = ind2sub(inds);
    }

    void sum(int dim){

        arma::umat aSubs = subs;
        aSubs.shed_row(dim);
        D.shed_row(dim);
        arma::uvec aInds = sub2ind(aSubs);

        inds = arma::unique(aInds);

        arma::Mat<typename T::elem_type> sumDataStore(inds.n_elem, 1, arma::fill::zeros);

        for (arma::uword c = 0; c < inds.n_elem; c++) {
            arma::uvec idx = aInds == inds(c);
            sumDataStore(c) = arma::sum(dataStore.elem(arma::find(idx)));
        }

        dataStore = sumDataStore;
        D.insert_rows(dim, arma::ones<arma::urowvec>(1));
        subs = ind2sub(inds);

        //We should always be sorted, but just in case...
        if (!inds.is_sorted()) {
            arma::uvec sortI = arma::sort_index(inds);
            inds = inds.elem(sortI);
            subs = subs.cols(sortI);
            dataStore = dataStore.rows(sortI);
        }
    }

    void squeeze(){
        D.elem(arma::find(D != 1));
        subs = ind2sub(inds);
    }

    arma::mat real(){ return arma::real(dataStore); };
    arma::mat imag(){ return arma::imag(dataStore); };

    arma::Col<typename T::elem_type> getAll() {
        return dataStore.elem(inds);
    }

    arma::uvec getIndex() {
        return inds;
    }

    arma::umat getSubs() {
        return subs;
    }

    arma::umat ind2sub(arma::uvec inSubs) {
        arma::umat retMat(D.n_elem, inSubs.n_elem);
        arma::uvec k = arma::cumprod(D);
        for (arma::uword i = D.n_elem -1; i > 0; i--){
            arma::uvec a = inSubs;
            arma::uword b = k(i -1);
            arma::uvec vi = a - b* arma::floor(a / b); //Remainder after division.
            arma::uvec vj = (inSubs - vi) / k(i - 1);
            retMat.row(i) = vj.t();
            inSubs = vi;
        }
        retMat.row(0) = inSubs.t();
        return retMat;
    }

    arma::uvec sub2ind(arma::umat inInd) {
        arma::urowvec retVec(inInd.n_cols, arma::fill::zeros);
        arma::uvec k = arma::cumprod(D);
        retVec += inInd.row(0);
        retVec += inInd.row(1)*D(0);
        for (arma::uword i = 2; i < D.n_elem; i++){
            retVec += k(i - 1) * inInd.row(i);
        }
        return retVec.t();
    }
};


class ndmat: public ndMat<arma::mat> {
public:
    ndmat(double *inVec, int *dims, int nDims) {
        arma::uvec thisD = arma::conv_to<arma::uvec>::from(arma::Col<int>(&(dims[0]), nDims));
        ndmat_base(arma::mat(&(inVec[0]), arma::prod(thisD), 1, false), thisD);
    };

    ndmat(arma::mat inMat, int *dims, int nDims) {
        ndmat_base(inMat, arma::conv_to<arma::uvec>::from(arma::Col<int>(&(dims[0]), nDims)));
    };

private:
    void ndmat_base(arma::mat inMat, arma::uvec thisD){
        D = thisD;
        dataStore = inMat;
        inds = arma::linspace<arma::uvec>(0, dataStore.n_elem - 1, dataStore.n_elem);
        subs = ind2sub(inds);
    }
};

class ndcx_mat: public ndMat<arma::cx_mat> {
public:
    ndcx_mat(double *inVecReal, double *inVecImag, int *dims, int nDims) {

        arma::uvec thisD = arma::conv_to<arma::uvec>::from(arma::Col<int>(&(dims[0]), nDims));
        ndcx_base(arma::cx_mat(
                arma::vec(&(inVecReal[0]), arma::prod(thisD), false),
                arma::vec(&(inVecImag[0]), arma::prod(thisD), false)),thisD);
    };

    ndcx_mat(arma::vec inVecReal, arma::vec inVecImag, int *dims, int nDims) {
        arma::uvec thisD = arma::conv_to<arma::uvec>::from(arma::Col<int>(&(dims[0]), nDims));
        ndcx_base(arma::cx_mat(inVecReal, inVecImag), thisD);
    };

    ndcx_mat(arma::cx_mat inMat, int *dims, int nDims) {
        ndcx_base(inMat,arma::conv_to<arma::uvec>::from(arma::Col<int>(&(dims[0]), nDims)));
    };

private:
    void ndcx_base(arma::cx_mat inMat, arma::uvec thisD){
        D = thisD;
        dataStore = inMat;
        inds = arma::linspace<arma::uvec>(0, dataStore.n_elem - 1, dataStore.n_elem);
        subs = ind2sub(inds);
    }
};