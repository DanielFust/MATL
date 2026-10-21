#ifndef _MATL_MATRIX_TRANSPOSE_PROTOTYPE__HPP
#define _MATL_MATRIX_TRANSPOSE_PROTOTYPE__HPP

namespace matl 
{
   /***********************************************************************************************
    * class _MATRIX_TRANSPOSE_
    * 
    * Class prototype for matrix transpose expressions.
    *************************************************************************************************/
    template<typename MAO>
    class _MATRIX_TRANSPOSE_ : public _EXPRESSION_<_MATRIX_TRANSPOSE_<MAO> >,
                               public _ELEMENTWISE_<_MATRIX_TRANSPOSE_<MAO> >
    {
      private:
        typename member_type<MAO>::elementwise_expression mat_;
      
      public:  
        // Resultant matrix data type (e.g. double, int, float, etc.)
        typedef typename promote<typename type_traits<MAO>::dtype>::type dtype;
 
        // Constructor (elementwise object)
        template<typename MAO_               = MAO,
                 enable_if_elementwise<MAO_> = true>
        _MATRIX_TRANSPOSE_(MAO_ const&);

        // Constructor (non-elementwise object)
        template<typename MAO_                = MAO,
                 disable_if_elementwise<MAO_> = true>
        _MATRIX_TRANSPOSE_(MAO_ const&);

        // Get number of rows of resultant matrix
        index_t rows() const;

        // Get number of columns of resultant matrix
        index_t cols() const;

        // Get number of elements in resultant matrix
        index_t size() const;

        // Get the element at the i-th row and j-th column
        dtype operator() (index_t i, index_t j) const;

        // Get the element at the i-th index
        dtype operator() (index_t i) const;

        // Check if evaluation target is safe
        template<typename MAT>
        constexpr bool unsafe_evaluation_target(_MATRIX_<MAT> const&) const;

        // Check if expression references address of target matrix
        template<typename MAT>
        constexpr bool references_address(_MATRIX_<MAT> const*) const;

        // Evaluate the expression to the target matrix
        template<typename MAT>
        void eval(_MATRIX_<MAT>&) const;
    };
}

#endif