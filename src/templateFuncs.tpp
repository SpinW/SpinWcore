//
// Created by ward_s on 15/03/18.
//

#include <armadillo>


template<typename T>
arma::Mat<typename T::elem_type> as_Mat(const T& X)
{
    return {X};
}


template<typename T>
T armaMod(const T &X, int n)
{
    const auto& an = as_Mat(X/n);
    return X - floor(an)*n;
}

template<typename T>
static arma::Cube<typename T::elem_type> permute(const T &cube, std::tuple<int, int, int> &order) {
    using namespace arma;

    int op1 = std::get<0>(order);
    int op2 = std::get<1>(order);
    int op3 = std::get<2>(order);

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
        if (perm == 213) {
            output.slice(s) = this_slice.t(); // cxr = (rxc)'
        } else {
            for (int r = 0; r < cube.n_rows; ++r) {
                if (perm == 231) {
                    output.slice(r)(arma::span::all, s) = vectorise(this_slice(r, arma::span::all));
                } else if (perm == 321) {
                    output.slice(r)(s, arma::span::all) = vectorise(this_slice(r, arma::span::all));
                } else {
                    for (int c = 0; c < cube.n_cols; ++c) {
                        switch (perm) {
                            case 132: {
                                output(r, s, c) = cube(r, c, s);
                            }
                                break;
                            case 312: {
                                output(s, r, c) = cube(r, c, s);
                            }
                                break;
                        }
                    }
                }
            }
        }
    }

    return output;
}