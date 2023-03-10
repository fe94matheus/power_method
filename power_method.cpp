#include "matrix.hpp"
#include "vector.hpp"
#include <iostream>
#include <tuple>
#include <cmath>

namespace flib
{
template<class T>//matrix<T> pwr(matrix<T> &a, vector<T> v)
std::tuple<vector<T>, T> pwr(matrix<T> &a, vector<T> v)
{
    matrix<T> z(a.rows(), 1);
    matrix<T> q0(a.rows(), 1);
    matrix<T> q0t(1, a.rows());
    vector<T> r(a.rows());


    T lambda0, lambda;

    q0 = v;
    q0t = v;

    double tol = 10e-16;

    do
    {
        lambda0 = (q0t * a * q0)[0];
        z = a * q0;
        q0 = z * (1.0/(z.two_norm()));
        for (uint64_t i = 0; i < q0.rows(); i++)
        {
            q0t[i] = q0[i];
        }
        lambda = (q0t * a * q0)[0];
        
    } while (std::abs(lambda - lambda0) > tol);
    r = q0;

    return std::make_tuple(r, lambda);
}
} // namespace flib
