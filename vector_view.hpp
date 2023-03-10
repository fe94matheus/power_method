#pragma once
#include <iostream>
#include <cmath>
#include <cassert>

namespace flib
{
template< class T>
class vector;

template< class T >
class vector_view
{
    private:
        T* data_;
        uint64_t size_;
        int64_t stride_;
    public:
        constexpr vector_view(){}
        
        constexpr vector_view( T* data, uint64_t size, int64_t stride )
        :data_( data ), size_( size ), stride_( stride ){}

        T &at( uint64_t i ) const
        {
            assert( i < size_);
            return data_[ i * stride_ ];
        }

        T &operator[]( uint64_t i ) const
        {
            return at( i );
        }
        
        constexpr uint64_t size() const
        {
            return size_;
        }

        T two_norm() const
        {
            T value;
            for(uint64_t i = 0; i < size(); ++i)
                value += at(i) * at(i);

            return sqrt( value );
        }

        friend std::ostream &operator<<( std::ostream &os, vector_view<T> const& v )
        {
            os << std::showpos << std::scientific;
            for (uint64_t i = 0; i < v.size(); ++i)
            {
                os << v.at(i) << " ";
            }
            os << std::endl;

            return os;
            
        }

        vector< T > operator+( vector_view< T > v ) const
        {
            assert( size_ == v.size_ );

            vector< T > u( size_ );

            for( uint64_t i = 0; i < size_; ++i )
                u[ i ] = data_[ i ] + v[ i ];

            return u;
        }

        vector_view< T >& operator=( vector_view< T > v )
        {
            
            assert( size_ == v.size() );
            for (uint64_t i = 0; i < size_; i++)
            {
                at(i) = v.at(i);
            }
            
            return *this;
        }

};
} // namespace flib
