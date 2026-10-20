/**************************************************************************************************
 * File: _matl_matrix_plus_matrix_implementation_.hpp
 * 
 * Author: Daniel Fust
 * Date: 11/1/2024
 * 
 * Class implementation for generic matrix-matrix sum expression template
 **************************************************************************************************/
namespace matl
{
   /**********************************************************************************************
    * Constructor
    * 
    * Initialize matrix-matrix sum expression from two _ELEMENTWISE_ objects
    * 
    * Parameters:
    *     mat1         Matrix or elementwise expression 1
    *     mat2         Matrix or elementwise expression 2
    **********************************************************************************************/
    template<typename MAO1, typename MAO2>
    template<typename MAO_1, typename MAO_2,
             enable_if_elementwise<MAO_1>, 
             enable_if_elementwise<MAO_2> >
    _MATRIX_PLUS_MATRIX_<MAO1,MAO2>::_MATRIX_PLUS_MATRIX_(MAO_1 const& mat1, 
                                                          MAO_2 const& mat2) : 
        mat1_(unqualified_cast(mat1)), 
        mat2_(unqualified_cast(mat2))
    {}

   /**********************************************************************************************
    * Constructor (elementwise - non-elementwise)
    * 
    * Initialize matrix-matrix sum expression from a matrix or elementwise expression and a 
    * non-elementwise expression.
    * 
    * Parameters:
    *    mat1         Matrix or elementwise expression
    *    expr2        Non-elementwise expression
    **********************************************************************************************/
    template<typename MAO1, typename MAO2>
    template<typename MAO_1, typename MAO_2,
             enable_if_elementwise<MAO_1>, 
             disable_if_elementwise<MAO_2> >
    _MATRIX_PLUS_MATRIX_<MAO1,MAO2>::_MATRIX_PLUS_MATRIX_(MAO_1 const& mat1, 
                                                          MAO_2 const& expr2) : 
        mat1_(unqualified_cast(mat1))
    {
        MemoryManager::allocate(this->mat2_, expr2.rows(), expr2.cols());
        expr2.eval(this->mat2_);
    }

   /**********************************************************************************************
    * Constructor (non-elementwise - elementwise)
    * 
    * Initialize matrix-matrix sum expression from a non-elementwise expression and an elementwise
    * expression or matrix.
    * 
    * Parameters:
    *     expr1       Non-elementwise expression
    *     mat2        Elementwise expression or matrix
    ***********************************************************************************************/
    template<typename MAO1, typename MAO2>
    template<typename MAO_1, typename MAO_2,
             disable_if_elementwise<MAO_1>, 
             enable_if_elementwise<MAO_2> >
    _MATRIX_PLUS_MATRIX_<MAO1,MAO2>::_MATRIX_PLUS_MATRIX_(MAO_1 const& expr1, 
                                                          MAO_2 const& mat2) : 
        mat2_(unqualified_cast(mat2))
    {
        MemoryManager::allocate(this->mat1_, expr1.rows(), expr1.cols());
        expr1.eval(this->mat1_);
    }  

   /**********************************************************************************************
    * Constructor (non-elementwise - non-elementwise)
    * 
    * Initialize matrix-matrix sum expression from two non-elementwise expressions.
    * 
    * Parameters:
    *     expr1       Non-elementwise expression 
    *     expr2       Non-elementwise expression
    ***********************************************************************************************/
    template<typename MAO1, typename MAO2>
    template<typename MAO_1, typename MAO_2,
             disable_if_elementwise<MAO_1>, 
             disable_if_elementwise<MAO_2> >
    _MATRIX_PLUS_MATRIX_<MAO1,MAO2>::_MATRIX_PLUS_MATRIX_(MAO_1 const& expr1, 
                                                          MAO_2 const& expr2)
    {   
        MemoryManager::allocate(this->mat1_, expr1.rows(), expr1.cols());
        MemoryManager::allocate(this->mat2_, expr2.rows(), expr2.cols());
        expr1.eval(this->mat1_);
        expr2.eval(this->mat2_);
    } 

   /**********************************************************************************************
    * Get indexed result of expression
    * 
    * Retrieve the resultant value at specified index pair
    * 
    * Template Parameters:
    *     MAO1          Matrix expression type 1
    *     MAO2          Matrix expression type 2
    * 
    * Parameters:
    *     i                Row index
    *     j                Column index 
    * Returns:
    *     dtype            Value of resultant matrix at index pair  
    **********************************************************************************************/
    template<typename MAO1, typename MAO2>
    typename _MATRIX_PLUS_MATRIX_<MAO1,MAO2>::dtype 
    _MATRIX_PLUS_MATRIX_<MAO1,MAO2>::operator() (index_t i, index_t j) const
    {return this->mat1_(i,j) + this->mat2_(i,j);}

   /**********************************************************************************************
    * Get single-indexed result of expression
    * 
    * Retrieve the resultant value at specified index
    * 
    * Template Parameters:
    *     MAO1          Matrix expression type 1
    *     MAO2          Matrix expression type 2
    * 
    * Parameters:
    *     i                Index of the element to get
    * Returns:
    *     dtype            Value of resultant matrix at index
    **********************************************************************************************/
    template<typename MAO1, typename MAO2>
    typename _MATRIX_PLUS_MATRIX_<MAO1,MAO2>::dtype
    _MATRIX_PLUS_MATRIX_<MAO1,MAO2>::operator() (index_t i) const
    {return this->mat1_(i) + this->mat2_(i);} 

   /**********************************************************************************************
    * Get number of rows of resultant matrix
    * 
    * Number of rows in expression. Downcasts to derived type.
    * 
    * Template Parameters:
    *     MAO1          Matrix expression type 1
    *     MAO2          Matrix expression type 2
    * 
    * Parameters:
    * Returns:
    *    index_t          Number of rows in resultant matrix
    **********************************************************************************************/
    template<typename MAO1, typename MAO2>
    index_t _MATRIX_PLUS_MATRIX_<MAO1,MAO2>::rows() const
    {return mat1_.rows();}

   /**********************************************************************************************
    * Get number of columns of resultant matrix
    * 
    * Number of columns in expression. Downcasts to derived type.
    * 
    * Template Parameters:
    *     MAO1          Matrix expression type 1
    *     MAO2          Matrix expression type 2
    * 
    * Parameters:
    * Returns:
    *    index_t          Number of columns in resultant matrix
    **********************************************************************************************/
   template<typename MAO1, typename MAO2>
   index_t _MATRIX_PLUS_MATRIX_<MAO1,MAO2>::cols() const
   {return mat1_.cols();}



   /**********************************************************************************************
    * Get number of elements in of resultant matrix
    * 
    * Number of elements in expression. Downcasts to derived type.
    * 
    * Template Parameters:
    *     MAO1          Matrix expression type 1
    *     MAO2          Matrix expression type 2
    * 
    * Parameters:
    * Returns:
    *    index_t          Number of elements in resultant matrix
    **********************************************************************************************/
    template<typename MAO1, typename MAO2>
    index_t _MATRIX_PLUS_MATRIX_<MAO1,MAO2>::size() const
    {return mat1_.size();}

   /***********************************************************************************************
    * bool unsafe_evaluation_target(_MATRIX_ const* target)
    *
    * Parameters:
    *     target       Candidate target matrix for evaluation of expression
    * Returns:
    *     'true' if the evaluation of the expression would corrupt memory.
    *     'false' otherwise.
    ************************************************************************************************/
    template<typename MAO1, typename MAO2>
    template<typename MAT>
    constexpr bool _MATRIX_PLUS_MATRIX_<MAO1,MAO2>::unsafe_evaluation_target(_MATRIX_<MAT> const& target) const
    {
        const bool branch1_unsafe = ((matl::references_address(&target,this->mat1_) && transposed_matrix(this->mat1_)) 
                                    || matl::unsafe_evaluation_target(target,this->mat1_));
        const bool branch2_unsafe = ((matl::references_address(&target,this->mat2_) && transposed_matrix(this->mat2_))
                                    || matl::unsafe_evaluation_target(target,this->mat2_));
        return (branch1_unsafe || branch2_unsafe);
    } 


   /***********************************************************************************************
    * bool references_address(_MATRIX_ const* target)
    *
    * Parameters:
    *     target       Candidate target matrix for evaluation of expression
    * Returns:
    *     'true' if the expression references the address of the target matrix.
    *     'false' otherwise.
    ************************************************************************************************/
    template<typename MAO1, typename MAO2>
    template<typename MAT>
    constexpr bool _MATRIX_PLUS_MATRIX_<MAO1,MAO2>::references_address(_MATRIX_<MAT> const* target) const
    {
        const bool branch1_ref = matl::references_address(target,this->mat1_);
        const bool branch2_ref = matl::references_address(target,this->mat2_);
        return (branch1_ref || branch2_ref);
    } 




   /************************************************************************************************
    * void eval(target) const
    * 
    * Evaluates the matrix sum expression to the target matrix
    * 
    * Template Parameters:
    *    MAO1          Matrix expression type 1
    *    MAO2          Matrix expression type 2
    *    MAT           Target matrix type
    * Parameters:
    *   target          Target matrix
    ************************************************************************************************/ 
    template<typename MAO1, typename MAO2>
    template<typename MAT>
    void _MATRIX_PLUS_MATRIX_<MAO1,MAO2>::eval(_MATRIX_<MAT>& target) const
    {
        if (is_same_storage_order<MAT,MAO1,MAO2>::value)
        {
            const index_t N = this->mat1_.size();
            for (index_t i=1; i<=N; ++i)
                target(i) = (*this)(i);
        }
        else
        {
            const index_t M = this->mat1_.rows();
            const index_t N = this->mat1_.cols();
            for (index_t i=1; i<=M; ++i)
                for (index_t j=1; j<=N; ++j)
                    target(i,j) = (*this)(i,j);
        }
 
        return;
    }
}