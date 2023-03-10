#pragma once
#include "vector_view.hpp"
#include "matrix_view.hpp"
#include "vector.hpp"
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <cassert>

namespace flib
{

template < class T >
class matrix_view;

template < class T >
class vector_view;

template < class T >
class vector;

template < class T >
class matrix
{
  private:

    T *data_;
    uint64_t rows_;
    uint64_t cols_;
    uint64_t alloc;

  public:

    ~matrix()
    {
      if( data_ )
        delete [] data_;
    }

    matrix(){}

    matrix( uint64_t rows, uint64_t cols )
    : data_( nullptr ), rows_( rows ), cols_( cols ), alloc( ( rows + 1 ) * ( cols + 1 ) )
    {
        data_ = new T[ alloc ];
    }

    matrix(matrix<T> &&a)
    : data_(a.data_), rows_(a.rows_), cols_(a.cols_)
    {
        a.data_ = nullptr;
    }

    uint64_t rows() const
    {
        return rows_;
    }
    uint64_t cols() const
    {
        return cols_;
    }

    T &at(uint64_t i) const
    {
        return data_[i];
    }

    T &operator[](uint64_t i) const
    {
        return at(i);
    }

    T &at(uint64_t rows, uint64_t cols) const
    {
        return data_[(cols_ * rows) + cols];
    }

    T &operator()(uint64_t rows, uint64_t cols)
    {
        return at(rows, cols);
    }

    void resize( uint64_t rows, uint64_t cols )
    {
        uint64_t new_alloc = rows * cols;
        if( new_alloc > alloc )
        {
            if( data_ )
            delete [] data_;
            alloc = std::max(new_alloc, 2 * alloc );
            data_ = new T[ alloc ];
        }
        rows_ = rows;
        cols_ = cols;
    }

    matrix<T> transpose()
    {
        matrix<T> t(cols_ , rows_);
        double* d = data_;

        for (int  i = 0; i < cols_; i++)
        {
            double* c = d;
            for (int j = 0; j< rows_; j++, c += cols_)
            {
                t(i,j) = *c;
            }
            d++;
        }

        return t;
    }

    void fill()
    {
        auto gen = []()
        {
            return std::rand() & 0xFFFF;
        };

        std::generate_n(data_, rows_ * cols_, gen);
    }

    matrix<T> eye()
    {

        assert( rows_ == cols_);
        matrix<double> I(rows_, rows_);

        for (int i = 0; i < rows_; ++i)
        {
            for (int j = 0; j < cols_; ++j)
            {
                I.data_[i*cols_ + j] = (i == j) ? 1.0 : 0.0;
            }
        }
        return I;
    }

    T two_norm()
    {
        T sum = 0;
        assert((rows_ == 1) || (cols_ == 1));
        for (uint64_t i = 0; i < rows_ * cols_; i++)
        {
            sum += data_[i] * data_[i];
        }
        return std::sqrt(sum);

    }

    bool issymmetric()
    {
        bool sign = true;
        for (uint64_t i = 0; i < rows_ * cols_; i++)
        {
            if(data_[i] != (*this).transpose()[i])
            {
                sign = false;
            }
        }
        return sign;
    }
    
    matrix<T> &operator+=(matrix<T> const &b)
    {
        assert((rows_ == b.rows()) && (cols_ == b.cols()));

        for (uint64_t k = 0; k < rows() * cols(); ++k)
            data_[k] += b.data_[k];

        return *this;
    }

    matrix<T> operator+(matrix<T> const &b) &&
    {
        return operator+=(b);
    }

    matrix<T> operator+(matrix<T> &b) const
    {
        assert((rows_ == b.rows()) && (cols_ == b.cols()));

        matrix<T> s(rows_, rows_);
        for (uint64_t k = 0; k < rows_ * cols_; ++k)
            s.data_[k] = data_[k] + b[k];
        return s;
    }

    matrix<T> operator*(T b) const
    {
        matrix<T> s(rows_, cols_);
        for (uint64_t k = 0; k < rows_ * cols_; ++k)
            s.data_[k] = data_[k] * b;
        return s;
    }

    matrix<T> operator-(matrix<T> b) const
    {
        assert((rows_ == b.rows()) && (cols_ == b.cols()));

        matrix<T> s(rows_, rows_);
        for (uint64_t k = 0; k < rows_ * cols_; ++k)
            s.data_[k] = data_[k] - b[k];
        return s;
    }
    matrix<T> operator+(matrix_view<T> b) const
    {
        assert((rows_ == b.rows()) && (cols_ == b.cols()));

        matrix<T> s(rows_, rows_);
        for (uint64_t k = 0; k < rows_ * cols_; ++k)
            s.data_[k] = data_[k] + b[k];
        return s;
    }

    matrix<T> operator*(matrix<T>& b)
    {
        assert( cols_ == b.rows() );

        matrix<T> m(rows_, b.cols());

        T* pm = m.data_;
        T const* pri;
        T const* pr = data_;
        for (uint64_t i = 0; i < rows_; ++i)
        {
            T const* pc = b.data_;
            for (uint64_t j = 0; j < b.cols(); ++j, ++pm, ++pc)
            {
                pri = pr;
                T const* pcj = pc;

                *pm = 0.0;
                for (int k = 0; k < cols_; ++k, ++pri, pcj += b.cols())
                {
                    *pm += *pri * *pcj;
                }
            }
            pr = pri;
        }
        return m;
    }

    matrix< T > operator-( matrix_view< T const > b ) const
    {
      assert((rows_ == b.rows()) && (cols_ == b.cols()));

      matrix< T > s( rows_, rows_ );

      for( uint64_t k = 0; k < rows_ * cols_; ++k )
          s.data_[k] = data_[k] - b[k];
      return s;
    }

    matrix< T > operator*( matrix_view<T> const& b ) const
    {
        assert( cols_ == b.rows() );

        matrix< T > m(rows_, b.cols());

        T* pm = m.data_;
        T const* pri;
        T const* pr = data_;
        for (uint64_t i = 0; i < rows_; ++i)
        {
            for (uint64_t j = 0; j < b.cols(); ++j, ++pm)
            {
                pri = pr;

                *pm = 0.0;
                for (int k = 0; k < cols_; ++k, ++pri)
                {
                    *pm += *pri * b.at(k,j);
                }
            }
            pr = pri;
        }
        return m;
    }

    /*
    matrix<T> &operator=( matrix_view<T const> mv)
    {
      resize( mv.rows(), mv.cols() );
      std::copy_n( mv.data(), rows_ * cols_, data_ );
      return *this;
    }

    matrix< T >& operator=( vector_view< T const>  v)
    {
        if( rows_ == 1 )
        {
            resize( 1, v.size() );
            std::copy_n( v.data(), cols_, data_);
        }
        else
        {
            assert( cols_ == 1 );
            resize( v.size(), 1 );
            std::copy_n( v.data(), rows_, data_);
        }

        return *this;
    }
    */
    matrix<T> &operator=(matrix<T> const& m)
    {
        assert( (rows_ == m.rows()) && (cols_ == m.cols()));
        for (uint64_t i = 0; i < rows_;  ++i)
        {
            for(uint64_t j = 0; j < cols_; ++j)
            {
                at(i,j) = m.at(i,j);
            }
        }

        return *this;
        
    }

    matrix<T> &operator=(matrix_view<T> const& mv)
    {
        assert( rows_ == mv.rows() && cols_ == mv.cols());
        for (uint64_t i = 0; i < rows_;  ++i)
        {
            for(uint64_t j = 0; j < cols_; ++j)
            {
                at(i,j) = mv.at(i,j);
            }
        }

        return *this;
        
    }

    matrix<T> &operator=(vector<T> const& v)
    {
        assert( (rows_ == v.size()) || (cols_ == v.size()));
        for (uint64_t i = 0; i < rows_ * cols_;  ++i)
        {
            at(i) = v.at(i);
        }

        return *this;
        
    }
    
    friend std::ostream& operator<<(std::ostream &os, matrix< T > const& m)
    {
        T const* a = m.data_;
        os << std::showpos << std::scientific;
        for (uint64_t i = 0; i < m.rows_; ++i)
        {
            for (uint64_t j = 0; j < m.cols_; ++j)
                os << *a++ << "  ";
            os << std::endl;
        }
        return os;
    }

    matrix_view< T > view(uint64_t i, uint64_t j, uint64_t rows, uint64_t cols, int32_t dr, uint32_t dc)
    {
        return matrix_view<T>( data_ + i * cols_ + j, rows, cols, dr * cols_ , dc );
    }

    matrix_view<T const> view(uint64_t i, uint64_t j, uint64_t rows, uint64_t cols, uint32_t dr, uint32_t dc) const
    {
        return matrix_view<T const>( data_ + i * cols_ + j, rows, cols, dr * cols_ , dc );
    }

    vector_view< T > col(uint64_t j, int begin = 0, int end = -1)
    {
        if (end < 0)
            end = rows_;

        return vector_view< T >(data_ + j + begin * cols_, end - begin, cols_);
    }
};


void test(uint64_t rows, uint64_t cols);
} // namespace felipe
