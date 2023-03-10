#pragma once

namespace flib
{

template < class T >
class matrix;

template < class T >
class matrix_view
{
  private:

    T *data_ = nullptr;
    uint64_t rows_ = 0;
    uint64_t cols_ = 0;
    uint64_t size_ = 0;
    int64_t row_stride_ = 0;
    int64_t col_stride_ = 0;

  public:

    matrix_view(){}

    matrix_view(T* data, uint64_t rows, uint64_t cols, int64_t row_stride = 1,
                int64_t col_stride = 1) : data_( data ), rows_( rows ), cols_( cols ),
                row_stride_( row_stride ), col_stride_( col_stride ){}


    constexpr uint64_t rows() const
    {
        return rows_;
    }

    constexpr uint64_t cols() const
    {
        return cols_;
    }

    T& at( uint64_t i ) const
    {
        return data_[ i ];
    }

    T& operator[](uint64_t i)
    {
        return at(i);
    }

    T& at( uint64_t row, uint64_t col ) const
    {
        return data_[ (row_stride_ * row ) + col * col_stride_ ];
    }

    T& operator()( uint64_t row, uint64_t col )
    {
        return at(row, col);
    }

    matrix_view< T >& operator=( matrix_view< T const > m )
    {
        assert( ( rows_ == m.rows() ) && ( cols_ == m.cols() ) );
        for (uint64_t i = 0; i < rows_; ++i)
        {
            for (uint64_t j = 0; j < cols_; ++j)
            {
                at(i, j) = m.at(i, j);
            }
            
        }
        return *this;
    }

    matrix_view<T> &operator=(matrix<T> const& m) 
    {
        assert( rows_ == m.rows());
        for (uint64_t i = 0; i < rows_; ++i)
        {
            for (uint64_t j = 0; j < cols_; ++j)
            {
                at(i, j) = m.at(i, j);
            }
            
        }

        return *this;
        
    }

    friend std::ostream& operator<<(std::ostream &os, matrix_view const& mv )
    {
      os << std::showpos;
      T const* a = mv.data_;
      for (uint64_t i = 0; i < mv.rows_; ++i, a += mv.row_stride_)
      {
          T const* r = a;
          for (uint64_t j = 0; j < mv.cols_; ++j, r += mv.col_stride_)
              os << *r << "  ";
          os << std::endl;
      }

      return os;
    }
};

} // namespace flib
