/*************************************************************************************************
 * _matl_matrix_matrix_sum_prototype_.hpp
 * 
 * Author: Daniel Fust
 * Date: 11/1/2024
 * 
 * Class prototype for generic matrix-matrix sum expression template
 *************************************************************************************************/
namespace matl
{
   /**********************************************************************************************
    * class _MATRIX_PLUS_MATRIX_
    * 
    * Generic matrix-matrix sum expression template
    * 
    * Template Parameters:
    *     MAO1: Matrix or expression type
    *     MAO2: Matrix or expression type
    **********************************************************************************************/
    template<typename MAO1, typename MAO2>
    class _MATRIX_PLUS_MATRIX_ : public _EXPRESSION_<_MATRIX_PLUS_MATRIX_<MAO1,MAO2> >
    {
      private:
        typename member_type<MAO1>::elementwise_expression mat1_;
        typename member_type<MAO2>::elementwise_expression mat2_;
      public:  
        // Resultant matrix data type (e.g. double, int, float, etc.)
        typedef typename promote<typename type_traits<MAO1>::dtype,
                                 typename type_traits<MAO2>::dtype>::type dtype;

        // Constructor (elementwise + elementwise)
        template<typename MAO_1 = MAO1,
                 typename MAO_2 = MAO2,
                 enable_if_elementwise<MAO_1> = true,
                 enable_if_elementwise<MAO_2> = true >
        _MATRIX_PLUS_MATRIX_(MAO_1 const&, MAO_2 const&);

        // Constructor (elementwise + non-elementwise expression)
        template<typename MAO_1 = MAO1,
                 typename MAO_2 = MAO2,
                 enable_if_elementwise<MAO_1>   = true,
                 disable_if_elementwise<MAO_2>  = true >
        _MATRIX_PLUS_MATRIX_(MAO_1 const&, MAO_2 const&);

        // Constructor (non-elementwise + elementwise)
        template<typename MAO_1 = MAO1,
                 typename MAO_2 = MAO2,
                 disable_if_elementwise<MAO_1> = true,
                 enable_if_elementwise<MAO_2>  = true >
        _MATRIX_PLUS_MATRIX_(MAO_1 const&, MAO_2 const&);

        // Constructor (non-elementwise + non-elementwise)
        template<typename MAO_1 = MAO1,
                 typename MAO_2 = MAO2,
                 disable_if_elementwise<MAO_1> = true,
                 disable_if_elementwise<MAO_2> = true >
        _MATRIX_PLUS_MATRIX_(MAO_1 const&, MAO_2 const&);

        // Get indexed result of expression
        dtype operator() (index_t, index_t) const;

        // Get single-indexed result of expression
        dtype operator() (index_t) const;

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

        // Evaluate matrix summation to target matrix
        template<typename MAT_T>
        void eval(_MATRIX_<MAT_T>& ) const;

    };
} 


