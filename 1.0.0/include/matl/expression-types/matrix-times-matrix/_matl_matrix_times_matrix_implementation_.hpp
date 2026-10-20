/**************************************************************************************************
 * File: _matl_matrix_times_matrix_implementation_.hpp
 * 
 * Description: Implementations for the matrix times matrix expression template.
 **************************************************************************************************/
namespace matl
{
   /***********************************************************************************************
    * MATRIX - MATRIX Constructor 
    * 
    * Constructs a matrix-matrix multiplication expression from two Matrix objects.
    * 
    * Parameters:
    *     mat1       Matrix 1
    *     mat2       Matrix 2
    ***********************************************************************************************/
   template<typename MAO1, typename MAO2>
    template<typename MAT1, typename MAT2> 
    _MATRIX_TIMES_MATRIX_<MAO1,MAO2>::_MATRIX_TIMES_MATRIX_(_MATRIX_<MAT1> const& mat1, 
                                                            _MATRIX_<MAT2> const& mat2) : 
        mat1_(unqualified_cast(mat1)), 
        mat2_(unqualified_cast(mat2)) 
    {
        #ifdef MATL_DEBUG
        std::cout << "_MATRIX_TIMES_MATRIX_: Matrix-Matrix Constructor" << std::endl;
        #endif
    }

   /***********************************************************************************************
    * MATRIX - EXPRESSION Constructor
    * 
    * Constructs a matrix-matrix multiplication expression from a Matrix object and an expression 
    * object.
    *  
    * Parameters:
    *     mat1        Matrix 
    *     expr2       Expression 
    ***********************************************************************************************/
    template<typename MAO1, typename MAO2>
    template<typename MAT1, typename EXPR2>
    _MATRIX_TIMES_MATRIX_<MAO1,MAO2>::_MATRIX_TIMES_MATRIX_(_MATRIX_<MAT1>      const& mat1, 
                                                            _EXPRESSION_<EXPR2> const& expr2) : 
        mat1_(unqualified_cast(mat1))                                                          
    {
        #ifdef MATL_DEBUG
        std::cout << "_MATRIX_TIMES_MATRIX_: Matrix-Expression Constructor" << std::endl;
        #endif
        MemoryManager::allocate(mat2_, expr2.rows(), expr2.cols());
        expr2.eval(mat2_);
    }


   /***********************************************************************************************
    * EXPRESSION - MATRIX Constructor
    *   
    * Constructs a matrix-matrix multiplication expression from an expression object and a Matrix
    * object.
    *
    * Parameters:
    *     expr1       Expression
    *     mat2        Matrix
    ***********************************************************************************************/
    template<typename MAO1, typename MAO2>
    template<typename EXPR1, typename MAT2>
    _MATRIX_TIMES_MATRIX_<MAO1,MAO2>::_MATRIX_TIMES_MATRIX_(_EXPRESSION_<EXPR1> const& expr1, 
                                                            _MATRIX_<MAT2>      const& mat2) : 
        mat2_(unqualified_cast(mat2))
    {
        #ifdef MATL_DEBUG
        std::cout << "_MATRIX_TIMES_MATRIX_: Expression-Matrix Constructor" << std::endl;
        #endif
        MemoryManager::allocate(mat1_, expr1.rows(), expr1.cols());
        expr1.eval(mat1_);
    } 

   /***********************************************************************************************
    * EXPRESSION - EXPRESSION Constructor
    * 
    * Constructs a matrix-matrix multiplication expression from two expression objects.
    * 
    * Parameters:
    *     expr1       Expression 1
    *     expr2       Expression 2
    ***********************************************************************************************/
   template<typename MAO1, typename MAO2>
    template<typename EXPR1, typename EXPR2>
    _MATRIX_TIMES_MATRIX_<MAO1,MAO2>::_MATRIX_TIMES_MATRIX_(_EXPRESSION_<EXPR1> const& expr1, 
                                                            _EXPRESSION_<EXPR2> const& expr2)
    {
        #ifdef MATL_DEBUG
        std::cout << "_MATRIX_TIMES_MATRIX_: Expression-Expression Constructor" << std::endl;
        #endif
        MemoryManager::allocate(mat1_, expr1.rows(), expr1.cols());
        MemoryManager::allocate(mat2_, expr2.rows(), expr2.cols());
        expr1.eval(mat1_);
        expr2.eval(mat2_);
    } 

   
    
   /***********************************************************************************************
    * index_t rows() const
    * 
    * Returns:
    *   index_t          Number of rows in resultant matrix
    ***********************************************************************************************/  
    template<typename MAO1, typename MAO2>        
    index_t _MATRIX_TIMES_MATRIX_<MAO1,MAO2>::rows() const
    {return mat1_.rows();}
     
   /***********************************************************************************************
    * index_t cols() const
    * 
    * Returns:
    *     index_t  Number of columns in resultant matrix
    ***********************************************************************************************/
    template<typename MAO1, typename MAO2>
    index_t _MATRIX_TIMES_MATRIX_<MAO1,MAO2>::cols() const
    {return mat2_.cols();}
    
   /***********************************************************************************************
    * index_t size() const
    * 
    * Returns:
    *     index_t          Number of elements in resultant matrix
    ***********************************************************************************************/
    template<typename MAO1, typename MAO2>
    index_t _MATRIX_TIMES_MATRIX_<MAO1,MAO2>::size() const
    {return mat1_.rows() * mat2_.cols();}

   /***********************************************************************************************
    * bool unsafe_evaluation_target(_MATRIX_ const* target)
    * 
    * Check if evaluation target is safe
    * 
    * Parameters:
    *     target           Candidate target matrix for evaluation of expression
    * Returns:
    *     'true' if the evaluation of the expression would corrupt memory.
    *     'false' otherwise.
    ***********************************************************************************************/
    template<typename MAO1, typename MAO2>
    template<typename MAT>
    constexpr bool _MATRIX_TIMES_MATRIX_<MAO1,MAO2>::unsafe_evaluation_target(_MATRIX_<MAT> const& target) const
    {return (matl::references_address(&target,this->mat1_) || matl::references_address(&target,this->mat2_));} 

   /***********************************************************************************************
    * bool references_address(_MATRIX_ const* target)
    * 
    * Check if the expression references the address of the target matrix.
    * 
    * Parameters:
    *     target           Matrix to check
    * Returns:
    *     'true' if the expression references the address of the target matrix
    *     'false' otherwise
    ***********************************************************************************************/
    template<typename MAO1, typename MAO2>
    template<typename MAT>
    constexpr bool _MATRIX_TIMES_MATRIX_<MAO1,MAO2>::references_address(_MATRIX_<MAT> const* target) const
    {return (matl::references_address(target,this->mat1_) || matl::references_address(target,this->mat2_));} 
    

   /***********************************************************************************************
    * void eval(_MATRIX_<MAT_T>& target) const
    * 
    * Evaluate the expression to the target matrix
    * 
    * Parameters:
    *     target           Target matrix
    ***********************************************************************************************/
    template<typename MAO1, typename MAO2>
    template<typename MAT>
    void _MATRIX_TIMES_MATRIX_<MAO1,MAO2>::eval(_MATRIX_<MAT>& target) const
    {
        const index_t M = mat1_.rows();
        const index_t N = mat2_.cols();
        const index_t K = mat1_.cols();

        for (index_t i = 1; i <= M; i++)
        {
            for (index_t j = 1; j <= N; j++)
            {
                typename type_traits<MAT>::dtype sum = 0;
                for (index_t k = 1; k <= K; k++)
                {
                    sum += mat1_(i,k) * mat2_(k,j);
                }
                target(i,j) = sum;
            }
        }
    }



   /***********************************************************************************************
    * member_type<MAO1>::matmul_expression const* mat1_address() const
    * 
    * This method may be used by _MARTIX_ objects during execution of the _MATRIX_ member 
    * 'operator=(_EXPRESSION_ const&)' to change execution behavior depending on whether the
    * lvalue matrix is also referenced in the expression.
    * 
    * Returns:
    *     member_type<MAO1>::matmul_expression*       Address of member matrix 1
    ***********************************************************************************************/
    template<typename MAO1, typename MAO2>
    constexpr void const* 
    _MATRIX_TIMES_MATRIX_<MAO1,MAO2>::mat1_address() const
    {return static_cast<void const*>(&(this->mat1_));}


   /***********************************************************************************************
    * member_type<MAO2>::matmul_expression const* mat2_address() const
    * 
    * This method may be used by _MARTIX_ objects during execution of the _MATRIX_ member 
    * 'operator=(_EXPRESSION_ const&)' to change execution behavior depending on whether the
    * lvalue matrix is also referenced in the expression.
    * 
    * Returns:
    *     member_type<MAO2>::matmul_expression*       Address of member matrix 2
    ***********************************************************************************************/
    template<typename MAO1, typename MAO2>
    constexpr void const*
    _MATRIX_TIMES_MATRIX_<MAO1,MAO2>::mat2_address() const
    {return static_cast<void const*>(&(this->mat2_));}

}