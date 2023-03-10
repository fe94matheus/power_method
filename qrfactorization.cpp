#include "vector.hpp"
#include "matrix.hpp"
#include <tuple>
#include <cmath>




namespace flib
{
template<class T>
std::tuple<vector<T>, T> house(vector_view<T> x)
{
    T beta;
    uint64_t n = x.size();
    
    vector<T> v(n);
    v = x;

    matrix<T> u(n, 1);
    matrix<T> ut(1, n);

    for (uint64_t j = 0; j < n;  ++j)
    {
        u[j] = v[j];
        ut[j] = u[j];
    }
    
    matrix<T> u1(n-1,1), ut1(1,n-1);

    u1 = u.view(1,0,n-1,1,1,1);
    ut1 = u.view(0,1,1,n-1,1,1);

    T sigma = (ut1 * u1)[0]; 
    
    v[0] = 1;

    if (sigma == 0)
        beta = 0;
    else
    {
        T mu = std::sqrt(x[0] * x[0] + sigma);
        if (x[0] <= 0)
            v[0] = x[0] - mu;
        else
        {
            v[0] = -sigma / (x[0] + mu);
        }
        beta = (2 * v[0] * v[0]) / (sigma + v[0] * v[0]);
        v = (1/v[0]) * v;
        
    }

    return std::make_tuple(v, beta); 
} 

template<class T>
vector<T> house(matrix<T> &a)
{
    vector<T> w(a.cols());
    T beta;

    for (uint64_t j = 0; j < a.cols();  ++j)
    {
        uint64_t m = a.rows()-j;
        uint64_t n = a.col(j,j).size(); 
        
        vector<T> v(n);
        matrix<T> I(m, m);

        std::tie(v, beta) = house(a.col(j,j));
        w[j] = beta;


        matrix<T> u(v.size(), 1);
        matrix<T> ut(1, v.size());

        u = v;
        ut = v;

        a.view( j, j, a.rows()-j, a.cols()-j, 1, 1) = ( I.eye() - ( u * ut ) * beta ) * 
        a.view( j, j, a.rows()-j, a.cols()-j, 1, 1);


        if(j < a.rows())
        {
            a.col(j, j+1) = v.view(1, v.size()-1);
        }
    }

    return w;
}

template<class T>
matrix<T> orthog_hhr(vector<T> const& v,  T beta)
{
    matrix<T> v1(v.size(), 1);
    matrix<T> vt1(1, v.size());
    matrix<T> I(v.size(), v.size());
    matrix<T> m(v.size(), v.size());

    v1 = v;
    vt1 = v;

    m = I.eye() - v1 * vt1 * beta;

    return m;
    
}

template<class T>
matrix<T> orthog_hhr(matrix<T> &a, vector<T> beta)
{
    matrix<T> q(a.rows(), a.rows());
    matrix<T> I(a.rows(), a.rows());
    matrix<T> q0(a.rows(), a.rows());
    
    vector<T> u(a.rows());

    q = I.eye();


    for (uint64_t j = 0; j < a.cols();  ++j)
    {
        u.assign(a.rows(), 0.0);
        u[j] = 1.0;
        for (uint64_t i = j+1; i < a.rows();  ++i)
        { 
            u[i] = a.col(j)[i];
        }        
        q0 = flib::orthog_hhr(u, beta[j]);

        q = q0 * q ;
        
    }
    
    return q;
}



} // namespace flib
