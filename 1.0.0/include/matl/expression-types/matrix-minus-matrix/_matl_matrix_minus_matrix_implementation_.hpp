/**************************************************************************************************
 * File: _matl_matrix_minus_matrix_implementation_.hpp
 * 
 * Description: Implementations for the matrix minus matrix expression template.
 **************************************************************************************************/
namespace matl
{
   /***********************************************************************************************
    * Constructor (elementwise - elementwise)
    * 
    * Initialize matrix-matrix subtraction expression from two elementwise expressions or matrices.
    * Arguments implement 'operator (index_t,index_t) const' interface so members 'mat1_' and 
    * 'mat2_' are initialized to the input expressions.
    * 
    * Member types:
    *     mat1_:    _ELEMENTWISE<EW> const&)
    *     mat2_:    _ELEMENTWISE<EW> const&)
    *   
    * 
    * Parameters:
    *     mat1       Matrix or expression 1
    *     mat2       Matrix or expression 2
    ***********************************************************************************************/
    template<typename MAO1, typename MAO2>
    template<typename MAO_1,
             typename MAO_2,
             enable_if_elementwise<MAO_1>, 
             enable_if_elementwise<MAO_2> >
    _MATRIX_MINUS_MATRIX_<MAO1,MAO2>::_MATRIX_MINUS_MATRIX_(MAO_1 const& mat1, 
                                                            MAO_2 const& mat2) : mat1_(unqualified_cast(mat1)), 
                                                                                 mat2_(unqualified_cast(mat2)) 
    {std::cout << "_MATRIX_MINUS_MATRIX_ elementwise - elementwise Constructor" << std::endl;} 


   /***********************************************************************************************
    * Constructor (elementwise - non-elementwise expression)
    * 
    * Initialize matrix-matrix subtraction expression from an elementwise expression or matrix and
    * a non-elementwise expression (e.g. matrix multiplication). 
    * 
    * Member types:
    *   mat1_:     _ELEMENTWISE_<EW> const&)
    *   mat2_:     _MATRIX_<MAT>
    * 
    * Parameters:
    *     mat1       Matrix or elementwise expression 
    *     expr2      non-elementwise expression 
    ***********************************************************************************************/
    template<typename MAO1, typename MAO2>
    template<typename MAO_1,
             typename MAO_2,
             enable_if_elementwise<MAO_1>, 
             disable_if_elementwise<MAO_2> >
    _MATRIX_MINUS_MATRIX_<MAO1,MAO2>::_MATRIX_MINUS_MATRIX_(MAO_1 const& mat1, 
                                                            MAO_2 const& expr2) : 
        mat1_(unqualified_cast(mat1)) 
    {
        std::cout << "_MATRIX_MINUS_MATRIX_ elementwise - non-elementwise Constructor" << std::endl;
        MemoryManager::allocate(this->mat2_,expr2.rows(), expr2.cols());
        expr2.eval(mat2_);

        // Call destructors on temporary objects
        const_cast<MAO_2 &>(expr2).~MAO_2();
    } 


   /***********************************************************************************************
    * Constructor (non-elementwise - elementwise)
    *
    * Initialize matrix-matrix subtraction expression from a non-elementwise expression and an
    * elementwise expression or matrix.
    * 
    * Member types:
    *     mat1_:     _MATRIX_<MAT>
    *     mat2_:     _ELEMENTWISE_<EW> const&)
    *
    * Parameters: 
    *     expr1      non-elementwise expression
    *     mat2       elementwise expression or matrix
    ***********************************************************************************************/
    template<typename MAO1, typename MAO2>
    template<typename MAO_1,
             typename MAO_2,
             disable_if_elementwise<MAO_1>, 
             enable_if_elementwise<MAO_2> >
    _MATRIX_MINUS_MATRIX_<MAO1,MAO2>::_MATRIX_MINUS_MATRIX_(MAO_1 const& expr1, 
                                                            MAO_2 const& mat2) : 
        mat2_(unqualified_cast(mat2))
    {
        std::cout << "_MATRIX_MINUS_MATRIX_ non-elementwise - elementwise Constructor" << std::endl;
        MemoryManager::allocate(this->mat1_, expr1.rows(), expr1.cols());
        expr1.eval(this->mat1_);

        // Call destructors on temporary objects
        const_cast<MAO_1 &>(expr1).~MAO_1();
    } 

   /***********************************************************************************************
    * Constructor (non-elementwise - non-elementwise)
    * 
    * Initialize matrix-matrix subtraction expression from two non-elementwise expressions.
    * 
    * Member types:
    *     mat1_:     _MATRIX_<MAT>
    *     mat2_:     _MATRIX_<MAT> 
    * 
    * Parameters:
    *     expr1       non-elementwise expression
    *     expr2       non-elementwise expression
    ***********************************************************************************************/
    template<typename MAO1, typename MAO2>
    template<typename MAO_1,
             typename MAO_2,
             disable_if_elementwise<MAO_1>, 
             disable_if_elementwise<MAO_2> >
    _MATRIX_MINUS_MATRIX_<MAO1,MAO2>::_MATRIX_MINUS_MATRIX_(MAO_1 const& expr1, 
                                                            MAO_2 const& expr2)
    {
        std::cout << "_MATRIX_MINUS_MATRIX_ non-elementwise - non-elementwise Constructor" << std::endl;
        MemoryManager::allocate(this->mat1_, expr1.rows(), expr1.cols());
        MemoryManager::allocate(this->mat2_, expr2.rows(), expr2.cols());
        expr1.eval(this->mat1_);
        expr2.eval(this->mat2_);

        // Call destructors on temporary objects
        const_cast<MAO_1 &>(expr1).~MAO_1();
        const_cast<MAO_2 &>(expr2).~MAO_2();
    } 
    

   /***********************************************************************************************
    * index_t rows() const
    * 
    * Get number of rows of resultant matrix
    * 
    * Number of rows in expression. Downcasts to derived type.
    * Returns:
    *   index_t          Number of rows in resultant matrix
    ***********************************************************************************************/  
    template<typename MAO1, typename MAO2>        
    index_t _MATRIX_MINUS_MATRIX_<MAO1,MAO2>::rows() const
    {return mat1_.rows();}
    
   /***********************************************************************************************
    * index_t cols() const
    * 
    * Get number of columns of resultant matrix
    * 
    * Number of columns in expression. Downcasts to derived type.
    * Returns:
    *     index_t  Number of columns in resultant matrix
    ***********************************************************************************************/
    template<typename MAO1, typename MAO2>
    index_t _MATRIX_MINUS_MATRIX_<MAO1,MAO2>::cols() const
    {return mat1_.cols();}


   /***********************************************************************************************
    * index_t size() const
    * 
    * Get number of elements in of resultant matrix
    * 
    * Number of elements in expression. Downcasts to derived type.
    * Returns:
    *     index_t          Number of elements in resultant matrix
    ***********************************************************************************************/
    template<typename MAO1, typename MAO2>
    index_t _MATRIX_MINUS_MATRIX_<MAO1,MAO2>::size() const
    {return mat1_.size();}


   /***********************************************************************************************
    * dtype operator() (index_t i, index_t j) const
    * 
    * Get indexed result of expression
    * 
    * Parameters:
    *     i           Row index
    *     j           Column index 
    * Returns:
    *    dtype       Value of resultant matrix at index pair
    ***********************************************************************************************/
    template<typename MAO1, typename MAO2>
    typename _MATRIX_MINUS_MATRIX_<MAO1,MAO2>::dtype
    _MATRIX_MINUS_MATRIX_<MAO1,MAO2>::operator() (index_t i, index_t j) const
    {return this->mat1_(i,j) - this->mat2_(i,j);} 

   /***********************************************************************************************
    * dtype operator() (index_t i) const
    * 
    * Get single-indexed result of expression
    * 
    * Parameters:
    *     i           Index of the element to get
    * Returns:
    *     dtype       Value of resultant matrix at index i 
    *************************************************************************************************/
    template<typename MAO1, typename MAO2>
    typename _MATRIX_MINUS_MATRIX_<MAO1,MAO2>::dtype
    _MATRIX_MINUS_MATRIX_<MAO1,MAO2>::operator() (index_t i) const
    {return this->mat1_(i) - this->mat2_(i);} 

   /***********************************************************************************************
    * bool unsafe_evaluation_target(_MATRIX_ const* target)
    * 
    * Check if evaluation of the expression to the target matrix would corrupt memory.
    * 
    * Parameters:
    *     target      Target matrix
    * Returns:
    *     'true' if the evaluation of the expression to target matrix would corrupt memory
    *     'false' otherwise
    ***********************************************************************************************/
    template<typename MAO1, typename MAO2>
    template<typename MAT>
    constexpr bool _MATRIX_MINUS_MATRIX_<MAO1,MAO2>::unsafe_evaluation_target(_MATRIX_<MAT> const& target) const
    {
        const bool branch1_unsafe = (matl::references_address(&target, this->mat1_) && transposed_matrix(this->mat1_))
                                    || matl::unsafe_evaluation_target(target, this->mat1_);
        const bool branch2_unsafe = (matl::references_address(&target, this->mat2_) && transposed_matrix(this->mat2_))
                                    || matl::unsafe_evaluation_target(target, this->mat2_);
        return (branch1_unsafe || branch2_unsafe);
    } 

   /***********************************************************************************************
    * bool references_address(_MATRIX_ const* target)
    *
    * Check if the expression references the address of the target matrix.
    *
    * Parameters:
    *     target      Target matrix
    * Returns:
    *     'true' if the expression references the address of the target matrix
    *     'false' otherwise
    ***********************************************************************************************/
    template<typename MAO1, typename MAO2>
    template<typename MAT>
    constexpr bool _MATRIX_MINUS_MATRIX_<MAO1,MAO2>::references_address(_MATRIX_<MAT> const* target) const
    {return (matl::references_address(target, this->mat1_) || matl::references_address(target, this->mat2_));} 

   /***********************************************************************************************
    * void eval(MAT& target) const
    * 
    * Evaluate matrix subtraction to target matrix
    * 
    * Parameters:
    *    target      Target matrix
    ***********************************************************************************************/
    template<typename MAO1, typename MAO2>
    template<typename MAT>
    void _MATRIX_MINUS_MATRIX_<MAO1,MAO2>::eval(_MATRIX_<MAT>& target) const
    {
        if (is_same_storage_order<MAT,MAO1,MAO2>::value)
        {
            const index_t N = this->mat1_.size();
            for (index_t i = 1; i <= N; i++)
                target(i) = (*this)(i);
        }
        else
        {
            const index_t M = this->mat1_.rows();
            const index_t N = this->mat1_.cols();
            for (index_t i = 1; i <= M; i++)
                for (index_t j = 1; j <= N; j++)
                    target(i,j) = (*this)(i,j);
        }
        
    } 
}