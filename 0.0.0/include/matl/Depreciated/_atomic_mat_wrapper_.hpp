#ifndef _MATL_ATOMIC_MAT_WRAPPER__HPP
#define _MATL_ATOMIC_MAT_WRAPPER__HPP
/**************************************************************************************************
 * Written by Daniel Fust
 *
 * This file contains the class definition and method implementations for the _ATOMIC_MAT_WRAPPER_
 * class which interfaces between a generic Matrix and _ATOMIC_EXPR_ (atomic expression) expreesion
 * templates.
 **************************************************************************************************/
#include "_matl_forward_declr_.hpp"
#include "_matl_type_traits_.hpp"

namespace matl
{    
   /***********************************************************************************************
    * class _ATOMIC_MAT_WRAPPER_
    * 
    * Row major atomic wrapper for a Matrix object. 
    * 
    * Template Arguments:
    *     MAT_T       Matrix type
    ***********************************************************************************************/
    template<typename MAT_T>
    class _ATOMIC_MAT_WRAPPER_ : public _ATOMIC_EXPR_< _ATOMIC_MAT_WRAPPER_<MAT_T> >
    {
      private:
        Matrix<MAT_T> const& mat_;

      public:
        // Matrix data type
        typedef typename TypeTraits<MAT_T>::type dtype;
          
        // Initialization Constructor
        _ATOMIC_MAT_WRAPPER_(Matrix<MAT_T> const& mat) : mat_(mat) {}

        // Destructor
        ~_ATOMIC_MAT_WRAPPER_() {}  

      protected:
        // Single index data getter
        dtype operator() (index_t i) const {return mat_(i);}

        // Double index data getter
        dtype operator() (index_t i, index_t j) const {return mat_(i,j);}

        // Evaluate expression
        //RowMatrix<dtype> eval() const;

        // Evaluate expression to target matrix
        template<typename MAT>
        void eval(Matrix<MAT>& target) const;

        // Number of rows in Matrix
        index_t rows() const {return mat_.rows();}

        // Number of columns in Matrix
        index_t cols() const {return mat_.cols();}

        // Number of elements in Matrix
        index_t size() const {return mat_.rows() * mat_.cols();}
            
    };



  


   /***********************************************************************************************
    * void _ATOMIC_MAT_WRAPPER_::eval(Matrix& target) const
    * 
    * "Evaluates" the matrix wrapper to a target matrix. This copies the contents of the Matrix 
    * wrapped by this _ATOMIC_MAT_WRAPPER_ object to a target matrix.
    * 
    * Template Parameters:
    *     MAT_T1 (implicit)    Matrix type of this wrapper
    *     MAT_T2 (implicit)    Matrix type of target
    * 
    * Parameters:
    *     target (Matrix&)     Target matrix for expression evaluation
    * 
    * Returns:
    ***********************************************************************************************/
    template<typename MAT_T1>
    template<typename MAT_T2>
    void _ATOMIC_MAT_WRAPPER_<MAT_T1>::eval(Matrix<MAT_T2>& target) const
    {   
        #ifdef DIM_CHECK // Check that target matrix is correctly sized
        if (  ((this->rows_ != target->rows()) || (this->rows_ != 0))
            ||((this->cols_ != target->cols()) || (this->cols_ != 0)) )
        {
            std::cerr << "Error in "
                      << "'void _ATOMIC_MAT_WRAPPER_<MAT_T1>::eval(Matrix<MAT_T2>* target)'"
                      << "\n    target Matrix with dimensions " << target->rows() << "x" 
                      << target->cols() << " incompatible with _ATOMIC_MAT_WRAPPER_ with "
                      << "dimensions " << this->rows() << "x" << this->cols() 
                      << std::endl;
            abort();          
        }    
        #endif

        // Evaluate expression to target matrix (CANDIDATE FOR CACHE-BLOCKING)
        for (index_t i = 0; i < this->rows(); ++i)
            for (index_t j = 0; j < this->cols(); ++j) 
                target(i,j) = (*this)(i,j);

        return;
    }
}    
#endif