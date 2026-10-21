/**************************************************************************************************
 * File: _matl_matrix_times_matrix_prototype_.hpp
 * 
 * Description: Prototype for the matrix times matrix expression template.
 **************************************************************************************************/
namespace matl
{
   /***********************************************************************************************
    * class _MATRIX_TIMES_MATRIX_
    * 
    * Generic matrix-matrix multiplication expression template. In order to reduce the number of
    * floating point operations, increase cache efficiency during evaluation, and provide 
    * compatibility with external BLAS implementations, expressions are evaluated and stored
    * as member matrices.
    * 
    * Template Parameters MAO* are of type 'Matrix' when MAO* is an expression (including
    * elementwise) and is of type 'const&' when MAO* is a _MATRIX_ type.
    ***********************************************************************************************/
    template<typename MAO1, typename MAO2>
    class _MATRIX_TIMES_MATRIX_ : public _EXPRESSION_<_MATRIX_TIMES_MATRIX_<MAO1,MAO2> >
    {
      private:
        typename member_type<MAO1>::matmul_expression mat1_;
        typename member_type<MAO2>::matmul_expression mat2_;
      public:  
        // Resultant matrix data type (e.g. double, int, float, etc.)
        typedef typename promote<typename type_traits<MAO1>::dtype,
                                 typename type_traits<MAO2>::dtype>::type dtype;

        // MATRIX - MATRIX Constructor 
        template<typename MAT1 = MAO1, typename MAT2 = MAO2>
        _MATRIX_TIMES_MATRIX_(_MATRIX_<MAT1> const&, _MATRIX_<MAT2> const&); 
        
        // MATRIX - EXPRESSION Constructor
        template<typename MAT1 = MAO1, typename EXPR2 = MAO2>
        _MATRIX_TIMES_MATRIX_(_MATRIX_<MAT1> const&, _EXPRESSION_<EXPR2> const&);

        // EXPRESSION - MATRIX Constructor
        template<typename EXPR1 = MAO1, typename MAT2 = MAO2>
        _MATRIX_TIMES_MATRIX_(_EXPRESSION_<EXPR1> const&, _MATRIX_<MAT2> const&);

        // EXPRESSION - EXPRESSION Constructor
        template<typename EXPR1 = MAO1, typename EXPR2 = MAO2>
        _MATRIX_TIMES_MATRIX_(_EXPRESSION_<EXPR1> const&, _EXPRESSION_<EXPR2> const&);

        // Get number of rows of resultant matrix
        index_t rows() const;

        // Get number of columns of resultant matrix
        index_t cols() const;

        // Get number of elements in resultant matrix
        index_t size() const;

        // Check if evaluation target is safe
        template<typename MAT>
        constexpr bool unsafe_evaluation_target(_MATRIX_<MAT> const&) const;

        // Check if expression references address of target matrix
        template<typename MAT>
        constexpr bool references_address(_MATRIX_<MAT> const*) const;

        // Evaluate matrix multiplication to target matrix
        template<typename MAT_T>
        void eval(_MATRIX_<MAT_T>& ) const;

        // Get address of member matrix 1
        //member_type<MAO1>::matmul_expression* mat1_address() const;
        constexpr void const* mat1_address() const;

        // Get address of member matrix 2
        //member_type<MAO2>::matmul_expression* mat2_address() const;
        constexpr void const* mat2_address() const;

    };

}