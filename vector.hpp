#pragma once
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <cassert>
#include "vector_view.hpp"
#include "matrix.hpp"

namespace flib
{

template<class T>
class vector_view;

template<class T>
class vector
{
    private:
        T* data_;
        uint64_t size_;
        uint64_t alloc_;
    public:
        /*~vector()
        {
            if( data_ )
                delete [] data_;
                std::cout << "Destructed!" << std::endl;
        }*/

        vector(){}

        vector(uint64_t size)
        :data_( nullptr ), size_( size ), alloc_( size )
        {
            data_ = new T[ size ];
        } 

        uint64_t size() const
        {
            return size_;
        }

        T& at( uint64_t i )
        {
            assert( i < size_) ;
            return data_[ i ];
        }

        T const& at( uint64_t i ) const
        {
            assert( i < size_) ;
            return data_[ i ];
        }

        T& operator[]( uint64_t i )
        {
            return at( i );
        }

        T& operator[]( uint64_t i ) const
        {
            return at( i );
        }

        void fill()
        {
            auto gen = []()
            {
                return std::rand() & 0xFFFF;
            };

            std::generate_n(data_, size_, gen);
        }

        template <class Gen>
        void generate( Gen gen )
        {
            std::generate_n( data_, size_, gen );
        }

        vector< T > append( vector_view< T const > v ) const
        {
            vector< T > a( size_ + v.size() );
            
            T *p = std::copy_n(data_, size_, a.data_);
            
            std::copy_n(v.data_, v.size_, p);
            
            return a;
        }

        void resize( uint64_t size )
        {
            if( size > alloc_ )
            {
                if( data_ )
                    delete [] data_;
                alloc_ = std::max( size , 2 * alloc_ );
                data_ = new T[ alloc_ ];
            }

            size_ = size;
        }

        void assign(uint64_t n, T value)
        {
            assert(n == size_);
            for (uint64_t i = 0; i < n;  ++i)
            {
                data_[i] = value;
            }
            
        }

        T two_norm() const
        {
            T value;
            for(uint64_t i = 0; i < size(); ++i)
                value += at(i) * at(i); 

            return sqrt( value );
        }

        operator vector_view< T >()
        {
            return vector_view< T >( data_, size_ );
        }

        operator vector_view< T const >() const
        {
            return vector_view< T const>( data_, size_ );
        }

        friend std::ostream &operator<<(std::ostream &os, vector<T> const& v)
        {
            T const* a = v.data_;
            os << std::showpos << std::scientific;
            for (uint64_t i = 0; i < v.size(); ++i)
            {
                os << *a++ << " ";
            }
            os << std::endl;

            return os;
            
        }

        vector< T > operator+( vector_view< T const > v ) const
        {
            assert( size_ == v.size() );
            vector<T> r(v.size());

            for( uint64_t i = 0; i < size_; ++i )
                r[i] = data_[i] + v[i];

            return r;
        }

        vector< T > operator-( vector_view< T const > v ) const
        {
            assert( size_ == v.size() );

            vector<T> r( v.size() );
            for (uint64_t i = 0; i < v.size(); ++i)
                r[i] = data_[i] - v[i];

            return r;
        }

        vector< T > operator*( T a ) const
        {
            vector< T > u( size_ );

            for( uint64_t i = 0; i < size_; i++ )
                u[ i ] = data_[ i ] * a;

            return u;
        }

        friend vector< T > operator*( T a, vector<T> const& v )  
        {
            vector< T > u( v.size_ );

            for( uint64_t i = 0; i < v.size_; i++ )
                u[ i ] = v.data_[ i ] * a;

            return u;
        }

        vector< T >& operator=( vector_view< T > const& v )
        {
            if( v.size() > alloc_ )
            {
                if( data_ )
                    delete [] data_;
                alloc_ = std::max(2 * alloc_, v.size() );
                data_ = new T[ alloc_ ];
            }

            size_ = v.size();
            for (uint64_t i = 0; i < size_; i++)
            {
                data_[i] = v[i];
            }
            return *this;
        }

        vector< T >& operator=( matrix< T > const& v )
        {
            if( v.rows() * v.cols() > alloc_ )
            {
                if( data_ )
                    delete [] data_;
                alloc_ = std::max(2 * alloc_, v.rows() * v.cols() );
                data_ = new T[ alloc_ ];
            }

            size_ = v.rows() * v.cols();
            for (uint64_t i = 0; i < size_; i++)
            {
                data_[i] = v[i];
            }
            return *this;
        }

        vector_view<T> view(uint64_t begin, int64_t end = -1, int64_t stride = 1)
        {
            return vector_view(data_ + begin, end + 1 - begin, stride);
        }

        vector_view< T const > view(uint64_t begin, int64_t end = -1, int64_t stride = 1) const
        {
            return vector_view< T const >(data_ + begin, end + 1 - begin, stride);
        }


};

void test(uint64_t n);
} // namespace flib
