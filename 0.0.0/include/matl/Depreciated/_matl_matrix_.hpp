#ifndef _MATL_MATRIX__HPP
#define _MATL_MATRIX__HPP

#include <iostream>
#include "_matl_forward_declr_.hpp"
#include "_matl_type_traits_.hpp"


namespace matl 
{

   /**********************************************************************************************
    * class Matrix
    * 
    * Base class for 'matl' matrix types.
    * 
    * Template Parameters:
    *     MAT_T       (Matrix) Matrix type (e.g. RowMatrix, ColMatrix, CSRMatrix, etc.)
    **********************************************************************************************/ 
    template<typename MAT_T>
    class Matrix
    {
      public:  

        // Matrix datatype
        typedef typename TypeTraits<MAT_T>::dtype dtype;

        /**************** Assignment & Copying *****************/
        template<typename MAT_T2>
        void operator= (Matrix<MAT_T2>& mat)
        {static_cast<MAT_T&>(*this)=mat;}

        /********************* Data Access *********************/

        // Single index getter
        dtype operator() (index_t i) const 
        {return static_cast<MAT_T&>(*this)(i);} 

        // Double index getter
        dtype operator() (index_t i, index_t j) const 
        {return static_cast<MAT_T&>(*this)(i,j);}

        // Single index accessor
        dtype& operator() (index_t i)
        {return static_cast<MAT_T&>(*this)(i);}

        // Double index accessor
        dtype& operator() (index_t i, index_t j)
        {return static_cast<MAT_T&>(*this)(i,j);}

        /********************* Utilities *********************/

        // Row Getter
        index_t rows() const
        {return static_cast<MAT_T&>(*this).rows();}

        // Coloumn Getter
        index_t cols() const
        {return static_cast<MAT_T&>(*this).cols();}

        // Size Getter
        index_t size() const
        {return static_cast<MAT_T&>(*this).size();}

        // Buffer size
        index_t buffer() const
        {return static_cast<MAT_T&>(*this).buffer();}
        

    };


#if 0


   /**********************************************************************************************
    * This class defines a Fortran-like slice mapping for Matrix operations
    **********************************************************************************************/ 
    class Slice
    {
        private:
            unsigned long start_;
            unsigned long stop_;
            std::size_t   stride_;
        public:
            
            /********************* Constructors *********************/
            
            // Default Constructor
            Slice();

            // Start-Stop slice with implicit stride of 1
            Slice(std::size_t start, std::size_t stop);

            // Fully defined constructor
            Slice(std::size_t start, std::size_t stop, std::size_t stride);

            /********************* Getters *********************/

            unsigned long start() const;

            unsigned long stop() const;

            unsigned long stride() const; 

            /********************* Setters *********************/

            void set_start(std::size_t) const;

            void set_stop(std::size_t) const;

            void set_stride(std::size_t) const;
    };




   /**********************************************************************************************
    * This is the foundational class for matrix expression templates
    **********************************************************************************************/
    template<typename E, class T>
    class _MATEXPRESSION_
    {
        protected:
            // Single Indexed Getter
            T operator() (unsigned long &i) const
            {return static_cast<E const&>(*this)(i);}

            // Double Indexed Getter
            T operator() (unsigned long &i, unsigned long &j) const
            {return static_cast<E const&>(*this)(i,j);}

            // Row Getter
            unsigned long rows()
            {return static_cast<E const&>(*this).rows();}

            // Column Getter
            unsigned long cols()
            {return static_cast<E const&>(*this).cols();}

            // Size Getter
            unsigned long size()
            {return static_cast<E const&>(*this).size();}

    };




   /**********************************************************************************************
    * Helper class template to map indices for sliced matrices
    **********************************************************************************************/
    //template<>





   


   /**********************************************************************************************
    * Matrix Summation expression template
    **********************************************************************************************/
    template<typename E1, typename E2, class T>
    class _MATSUM_ : public _MATEXPRESSION_<_MATSUM_<E1, E2, T>, T>
    {
        private:
            E1 const& _mat1;
            E2 const& _mat2;
        public:
            _MATSUM_(E1 const& mat1, E2 const& mat2) : _mat1{mat1}, _mat2{mat2}{}

            T operator() (unsigned long, unsigned long) const;

            T operator() (std::size_t, std::size_t) const;

            unsigned long size() const {return this->_mat1.size();}

            unsigned long rows() const {return this->_mat1.rows();}

            unsigned long cols() const {return this->_mat2.cols();}
    };












   /*=============================================================================================*
    *                        Slice:       CONSTRUCTORS
    *=============================================================================================*/
    
    
   /**********************************************************************************************
    * Default constructor
    *     encompases entire dimension of matrix with a stride of 1
    **********************************************************************************************/ 
    Slice::Slice() : start_{0}, stop_{std::numeric_limits<unsigned long>::max()}, stride_{1} {}



   /**********************************************************************************************
    * Contiguous subset constructor
    *     Start and stop indices with a stride of 1
    **********************************************************************************************/ 
    Slice::Slice(std::size_t start, std::size_t stop) : start_{start}, stop_{stop}, stride_{1} 
    {
        // Check for compatibility between unsigned index types and initializer inputs
        #ifdef UNSIGNED_COMP
        if (start < 0) || (stop > std::numeric_limits<unsigned long>::max())
        {
        -    std::cerr << "Error in 'Slice(std::size_t start, std::size_t stop)'"
                      << "\n    start=" << start << " and stop=" << stop 
                      << " values for Slice initialization are incompatible with matrix index " 
                      << "type 'unsigned long'" << std::endl;
            abort();          
        }
        #endif

        bool start_compatible = (start>0) && (start<=std::numeric_limits<unsigned long>::max());
        bool stop_compatible =  (stop<=std::numeric_limits<unsigned long>::max()) && (stop>0);
        
        this->start_ = start*(start_compatible) + (!start_compatible);
        this->stop_ =  stop*(stop_compatible) 
                     + std::numeric_limits<unsigned long>::max()*(!stop_compatible);
    }

    
   /**********************************************************************************************
    * Fully defined constructor
    **********************************************************************************************/ 
    Slice::Slice(std::size_t start, std::size_t stop, std::size_t stride) : stride_{stride}
    {
        //Check that the stride is not 0
        #ifdef MISC_DEBUG
        if (stride == 0)
        {
            std::cerr << "Error in 'Slice(std::size_t start, std::size_t stop, "
                      << "std::size_t stride)'"
                      << "\n    invalid stride=0"
                      << std::endl;
            abort();          
        }
        #endif
        // Check for compatibility between unsigned index types and initializer inputs
        #ifdef UNSIGNED_COMP
        if (start < 0) || (stop > std::numeric_limits<unsigned long>::max())
        {
        -    std::cerr << "Error in 'Slice(std::size_t start, std::size_t stop, "
                       << "std::size_t stride)'"
                       << "\n    start=" << start << " and stop=" << stop 
                       << " values for Slice initialization are incompatible with matrix index " 
                       << "type 'unsigned long'" << std::endl;
            abort();          
        }
        #endif

        bool start_compatible = (start>0) && (start<=std::numeric_limits<unsigned long>::max());
        bool stop_compatible =  (stop<=std::numeric_limits<unsigned long>::max()) && (stop>0);
        
        // Forward Stride
        if (stride > 0) 
        {
            this->start_ = start*(start_compatible) + 1*(!start_compatible);
            this->stop_ =  stop*(stop_compatible) 
                         + std::numeric_limits<unsigned long>::max()*(!stop_compatible);
        }
        // Backward Stride
        else
        {
            this->start_ =  start*(start_compatible) 
                          + std::numeric_limits<unsigned long>::max()*(!start_compatible);
            this->stop_ = stop*(stop_compatible) + (!stop_compatible);
        }
    }





   /**********************************************************************************************
    * Matrix Summation 
    **********************************************************************************************/
    template<typename E1, typename E2, class T>  
    inline T _MATSUM_<E1,E2,T>::operator() (unsigned long i, unsigned long j) const
    {
        return this->_mat1(i,j) + this->_mat2(i,j);
    }


   /**********************************************************************************************
    * Matrix Summation with index type preservation
    **********************************************************************************************/
    template<typename E1, typename E2, class T>  
    inline T _MATSUM_<E1,E2,T>::operator() (std::size_t i, std::size_t j) const
    {
        return this->_mat1(i,j) + this->_mat2(i,j);
    } 




   /**********************************************************************************************
    * Template Metaprogramming helper to determine arithmetic return type of two types
    **********************************************************************************************/
    template<class T1, class T2>
    constexpr decltype(auto) arithmetic_return_type()
    {
        T1 val1 = 0;
        T2 val2 = 0;
        return val1+val2;
    }
 #endif

}



#endif