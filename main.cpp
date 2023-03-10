#include "matrix.hpp"
#include "vector.hpp"
#include "qrfactorization.cpp"
#include "power_method.cpp"
#include <gtest/gtest.h>
#include <iostream>
#include <tuple>


TEST(QRFACTORIZATION, normCol)
{
    std::srand(1110);
    
    uint64_t m = 4;
    uint64_t n = 3;

    double tol = 10e-8;
    flib::matrix<double> a(m, n);
    flib::matrix<double> b(m, n);
    
    uint64_t samples = 1000000;
    for(uint64_t i = 0; i < samples; ++i)
    {
        a.fill();
        b = a;

        flib::house(a);

        for (uint64_t j = 0; j < n; ++j)
        {
            double value1, value2;

            value1 = b.col(j).two_norm();
            value2 = a.col(j, 0, j+1).two_norm();
            ASSERT_LE(std::abs (value1 - value2), tol);
        }
    }

    
}

TEST(QRFACTORIZATION, orthogHhr)
{
    std::srand(1110);
    
    uint64_t m = 4;
    uint64_t n = 3;

    double tol = 10e-8;
    
    uint64_t samples = 1000000;
    
    flib::matrix<double> a(m, n);
    flib::matrix<double> b(m, n);
    flib::matrix<double> q(m, m);
    flib::vector<double> beta(n);

    for (uint64_t i = 0; i < samples; ++i)
    {
        a.fill();
        b = a;
        beta = flib::house(a);
        q = flib::orthog_hhr(a, beta);

        for (uint64_t j = 0; j < a.cols(); ++j)
        {
            for (uint64_t k = j; k < a.cols(); ++k)
            {
                double value1, value2;
                value1 = a.at(j,k);
                value2 = (q * b).at(j,k);
                ASSERT_LE(std::abs(value1-value2),tol);
            } 
        }   
    }
}

/*TEST(POWERMETHOD, eigenvalue)
{
    uint64_t n = 2;
    flib::matrix<double> a(n, n);
    flib::vector<double> q(n);

    flib::vector<double> r(n);
    double eigen_value;

    a[0] = 5.0;
    a[1] = -2.0;
    a[2] = -2.0;
    a[3] = 8.0;

    a[0] = 2.0;
    a[1] = 1.0;
    a[2] = 1.0;
    a[3] = 2.0;
    
    q[0] = 1.0;
    q[1] = 1.0;

    if (a.issymmetric())
        std::cout << "It's symmetric" << std::endl;
    else
        std::cout << "It's not symmetric" << std::endl;

    
    std::tie(r, eigen_value) = flib::pwr(a, q);

    std::cout << "eigen_value: " << eigen_value << std::endl;
    std::cout << r;

}*/

TEST(POWERMETHOD, eigenvalue)
{
    uint64_t n = 2;
    
    flib::matrix<double> a(n, n);
    flib::matrix<double> x(a.rows(), 1);
    
    flib::vector<double> q(n);
    flib::vector<double> r(n);

    double eigen_value;
    
    double tol = 10e-8;
    uint64_t samples = 1;

    q.assign(n, 1.0);

    a[0] = 5.0;
    a[1] = -2.0;
    a[2] = -2.0;
    a[3] = 8.0;
    
    std::tie(r, eigen_value) = flib::pwr(a, q);
    
    x = r;

    for (uint64_t i = 0; i < n; i++)
    {
        ASSERT_LE(std::abs((a*x).at(i) - (x*eigen_value).at(i)), tol);
    }
    
}

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc,argv);
    return RUN_ALL_TESTS();

    return 0;
}