//
// Created by ward_s on 15/03/18.
//

#include "../include/templateFuncs.h"

template<typename T>
arma::Mat<typename T::elem_type> as_Mat(const T& X)
{
    return {X};
}


template<typename T>
T armaMod(const T &X, int n)
{
    const auto& a = as_Mat(X);
    return a - floor(a/n)*n;
}

template <typename T>
static arma::Cube<T> permute(arma::Cube<T>& cube, const std::tuple<arma::uword,arma::uword,arma::uword>& order)
{
    using namespace arma;
    uword idx1 = std::get<0>(order);
    uword idx2 = std::get<1>(order);
    uword idx3 = std::get<2>(order);

    uword perm = idx1*100 + idx2*10 + idx3;
    if (perm == 123) {
//        output(r, c, s)
        return cube;
    }


    u32_vec dimension({cube.n_rows, cube.n_cols, cube.n_slices});

    uword rows = dimension(idx1 - 1);
    uword cols = dimension(idx2 - 1);
    uword slis = dimension(idx3 - 1);

    Cube<T> output(rows, cols, slis, fill::zeros);

    for (arma::uword s = 0; s < cube.n_slices; ++s) {
        arma::mat sl = cube.slice(s);
        switch (perm) {
            case 132: {
//                output(r, s, c)
                output(arma::span::all, arma::span(s,s), arma::span::all) = sl;
            }
                break;
            case 213: {
//                output(c, r, s)
                output(arma::span::all, arma::span::all, arma::span(s,s)) = sl.t();
            }
                break;
            case 231: {
//                output(c, s, r)
                output(arma::span::all, arma::span(s,s), arma::span::all) = sl.t();
            }
                break;
            case 312: {
//                output(s, c, r)
                output(arma::span(s,s), arma::span::all, arma::span::all) = sl.t();
            }
                break;
            case 321: {
//                output(s, r c)
                output(arma::span(s,s), arma::span::all, arma::span::all) = sl;
            }
                break;
        }
    }

//    switch (perm)
//    {
//        case 123:
//        {
//            output = cube; // identity
//        }
//            break;
//        case 132:
//        {
//            for (int c = 0; c < cube.n_cols; ++c)
//                for (int r = 0; r < cube.n_rows; ++r)
//                    for (int s = 0; s < cube.n_slices; ++s)
//                        output(r, s, c) = cube(r, c, s);
//        }
//            break;
//        case 213:
//        {
//            for (int c = 0; c < cube.n_cols; ++c)
//                for (int r = 0; r < cube.n_rows; ++r)
//                    for (int s = 0; s < cube.n_slices; ++s)
//                        output(c, r, s) = cube(r, c, s);
//        }
//            break;
//        case 231:
//        {
//            for (int c = 0; c < cube.n_cols; ++c)
//                for (int r = 0; r < cube.n_rows; ++r)
//                    for (int s = 0; s < cube.n_slices; ++s)
//                        output(c, s, r) = cube(r, c, s);
//        }
//            break;
//        case 312:
//        {
//            for (int c = 0; c < cube.n_cols; ++c)
//                for (int r = 0; r < cube.n_rows; ++r)
//                    for (int s = 0; s < cube.n_slices; ++s)
//                        output(s, r, c) = cube(r, c, s);
//        }
//            break;
//        case 321:
//        {
//            for (int c = 0; c < cube.n_cols; ++c)
//                for (int r = 0; r < cube.n_rows; ++r)
//                    for (int s = 0; s < cube.n_slices; ++s)
//                        output(s, c, r) = cube(r, c, s);
//        }
//            break;
//    }

    return output;
}