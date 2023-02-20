#ifndef _BLAS_ROWMAT_ATOMIC__HPP
#define _BLAS_ROWMAT_ATOMIC__HPP
/**************************************************************************************************
 * Written by Daniel Fust
 *
 * This file contains the method declarations and implementations of the _ROWMAT_ATOMIC_ class
 **************************************************************************************************/
#include "_blas_forward_declr_.hpp"
#include "_blas_type_traits_.hpp"

namespace blas
{    
  /***********************************************************************************************
   * class _ROWMAT_ATOMIC_
   * 
   * Row major atomic wrapper for a RowMatrix object. 
   * 
   * Template Arguments:
   *     T          Matrix data type (e.g. double, int, etc.)
   ***********************************************************************************************/
   template<class T>
   class _ROWMAT_ATOMIC_ : public _ROWMAT_ATOMIC_EXPR_< _ROWMAT_ATOMIC_<T> >
   {
     private:
       RowMatrix<T> const& mat_;

     public:
       // Matrix data type
       typedef T dtype;
        
       // Initialization Constructor
       _ROWMAT_ATOMIC_(RowMatrix<T> const& mat) : mat_(mat) {}

       // Destructor
       ~_ROWMAT_ATOMIC_() {}  

       // Single index data getter
       dtype operator() (index_t i) const {return this->mat_(i);}

       // Double index data getter
       dtype operator() (index_t i, index_t j) const {return this->mat_(i,j);}

       // Evaluate expression
       RowMatrix<dtype> eval() const;

       // Evaluate expression to target matrix
       template<typename MAT_T>
       void eval(Matrix<MAT_T>& target) const;

       // Number of rows in Matrix
       index_t rows() const {return this->mat_.rows();}

       // Number of columns in Matrix
       index_t cols() const {return this->mat_.cols();}

       // Number of elements in Matrix
       index_t size() const {return this->mat_.rows() * this->mat_.cols();}
           
   };





   /***********************************************************************************************
    * void _ROWMAT_ATOMIC_::eval(Matrix& target) const
    * 
    * "Evaluates" the row matrix to a target matrix. This copies the contents of the RowMatrix 
    * wrapped by this _ROWMAT_ATOMIC_ object to a target matrix of any kind.
    * 
    * Template Parameters:
    *     T         (implicit) RowMatrix data type
    *     MAT_T     (implicit) Target matrix type 
    * 
    * Parameters:
    *     target    (Matrix*) target matrix
    * 
    * Returns:
    ***********************************************************************************************/
    template<class T>
    template<typename MAT_T>
    void _ROWMAT_ATOMIC_<T>::eval(Matrix<MAT_T>& target) const
    {
       
        #ifdef DIM_CHECK // Check that target matrix is correctly sized
        if (  ((this->rows_ != target->rows()) || (this->rows_ != 0))
            ||((this->cols_ != target->cols()) || (this->cols_ != 0)) )
        {
            std::cerr << "Error in 'void _ROWMAT_SUM_<E1,E2>::eval(Matrix<MAT_T>* target)"
                        "&expression)'"
                      << "\n    target Matrix with dimensions " << target->rows() << "x" 
                      << target->cols() << " incompatible with _ROWMAT_SUM_ with dimensions " 
                      << this->rows() << "x" << this->cols() 
                      << std::endl;
            abort();          
        }    
        #endif
        
        // Evaluate expression to target matrix
        for (index_t i = 0; i < this->rows(); ++i)
        {
            for (index_t j = 0; j < this->cols(); ++j) {
                target(i,j) = (*this)(i,j);
            }
        }

        return target;
    }
}    




#endif