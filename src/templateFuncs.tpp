//
// Created by ward_s on 15/03/18.
//

#include <armadillo>


template<typename T>
arma::Mat<typename T::elem_type> as_Mat(const T &X) {
    return {X};
}


template<typename T>
T armaMod(const T &X, int n) {
    const auto &an = as_Mat(X / n);
    return X - floor(an) * n;
}

template<typename T>
static arma::Cube<typename T::elem_type> permute(const T &cube, int order[]) {
    using namespace arma;

    int op1 = order[0];
    int op2 = order[1];
    int op3 = order[2];

    int perm = op1 * 100 + op2 * 10 + op3;

    if (perm == 123) {
        return cube;
    }

    uword dimension[] = {cube.n_rows, cube.n_cols, cube.n_slices};

    uword rows = dimension[op1 - 1];
    uword cols = dimension[op2 - 1];
    uword slis = dimension[op3 - 1];


    Cube<typename T::elem_type> output(rows, cols, slis, fill::zeros);

#pragma omp parallel for
    for (int s = 0; s < cube.n_slices; ++s) {
        auto this_slice = cube.slice(s); // rxc

        if (op1 == 2 & op2 == 1) {
            if (op3 == 3) { // 2 1 3
                output.slice(s) = this_slice.t(); // cxr = (rxc)'
            } else {
                throw std::runtime_error(std::string("Element is duplicated"));
            }
        } else { // [1 2 3] [1 2 3] [1 2 3]
            if (op3 == 1) { // [1 2 3] [1 2 3] [1]
                for (int r = 0; r < cube.n_rows; ++r) {
                    if (op1 == 2) { // [2] [1 2 3] [1]
                        output.slice(r)(arma::span::all, s) = vectorise(this_slice(r, arma::span::all));
                    } else if (op1 == 3) { // [3] [1 2 3] [1]
                        output.slice(r)(s, arma::span::all) = this_slice(r, arma::span::all);
                    } else {
                        throw std::runtime_error(std::string("Element is duplicated"));
                    }
                }
            } else if (op3 == 2) { // [1 2 3] [1 2 3] [2]
                for (int c = 0; c < cube.n_cols; ++c) {
                    if (op2 == 1) { // [1 2 3] [1] [2]
                        output.slice(c)(s, arma::span::all) = this_slice(arma::span::all, c).t();
                    } else if (op2 == 3) { // [1 2 3] [3] [2]
                        output.slice(c)(arma::span::all, s) = this_slice(arma::span::all, c);
                    } else {
                        throw std::runtime_error(std::string("Element is duplicated"));
                    }
                }
            } else{
                throw std::runtime_error(std::string("Element is duplicated"));
            }
        }
    }

    return output;
}

template<typename T>
arma::Mat<typename T::elem_type> mmat(const T &X, const T &Y, int dim[]) {

    arma::uword nA[2];
    arma::uword nB[2];

    if (dim[0]==2 && dim[1]==1){
        nA[0] = X.n_cols; nA[1] = X.n_rows;
        nB[0] = Y.n_cols; nB[1] = Y.n_rows;
    } else{
        nA[0] = X.n_rows; nA[1] = X.n_cols;
        nB[0] = Y.n_rows; nB[1] = Y.n_cols;
    }
    arma::uword repv[3] = {1, 1, 1};
    repv[dim[0]-1] = nA[0];

    int idx[3];
    idx[(dim[0]+1)%2] = 3;
    idx[dim[0]%2]     = dim[0];
    idx[2]            = dim[1];

    // Create temp and return type.
    arma::Cube<typename T::elem_type>A(X.n_rows, X.n_cols,nB[1], arma::fill::zeros);
    arma::Cube<typename T::elem_type>B;
    arma::Mat <typename T::elem_type>C;

    // Build A
    A.each_slice([X](arma::Mat<typename T::elem_type>&this_slice){this_slice = X;});

    // Build B
    if (idx[1] == 1)
    {
        B = arma::zeros(repv[0], Y.n_rows, Y.n_cols);
        for (arma::uword i = 0; i < B.n_slices; i++) {
            B.slice(i) = arma::repmat(Y(arma::span::all, i).t(), repv[0], repv[1]);
        }
    } else {
        B = arma::zeros(Y.n_cols, repv[1], Y.n_rows);
        for (arma::uword i = 0; i < B.n_slices; i++) {
            B.slice(i) = arma::repmat(Y(i,arma::span::all).t(), repv[0], repv[1]);
        }
    }

    // Inplace element multiplication A*B
    A%=B;

    // Create return matrix C
    if (dim[0] == 1) {
        // We sum along cols...
        C = arma::zeros(repv[0], A.n_slices);
        for (arma::uword i = 0; i < A.n_cols; i++){
            C += A(arma::span::all, arma::span(i,i), arma::span::all);
        }
    }
    else {
        // We sum along rows
        C = arma::zeros(A.n_slices, repv[1]);
        for (arma::uword i = 0; i < A.n_rows; i++){
            C += A(arma::span(i,i), arma::span::all, arma::span::all);
        }
    }
    return C;
}