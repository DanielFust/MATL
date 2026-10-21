#ifndef _MATL_MATRIX_DIVIDED_BY_SCALAR_PROTOTYPE__HPP
#define _MATL_MATRIX_DIVIDED_BY_SCALAR_PROTOTYPE__HPP

namespace matl 
{
   /***********************************************************************************************
    * class _MATRIX_DIVIDED_BY_SCALAR_
    * 
    * Expression template class for the matrix divided by scalar operation. This class is defined
    * in addition to the _MATRIX_TIMES_SCALAR_ class to preserve expected behavior of finite
    * precision arithmetic with types such as integers which have zero reciprocals.
    ***********************************************************************************************/ 
    template<typename MAO, class T>
    class _MATRIX_DIVIDED_BY_SCALAR_ : public _EXPRESSION_<_MATRIX_DIVIDED_BY_SCALAR_<MAO,T> >,
                                       public _ELEMENTWISE_<_MATRIX_DIVIDED_BY_SCALAR_<MAO,T> >
    {
      private:
        // Expression or matrix
        typename member_type<MAO>::elementwise_expression mat_;

        // Scalar (type promotion NOT applied to preserve expected behavior)
        T scalar_;  

      public:
        // Typedefs
        typedef typename promote<typename type_traits<MAO>::dtype, T>::type dtype;

        // Constructor (elementwise / scalar)
        template<typename MAO_               = MAO,
                 class    T_                 = T,
                 enable_if_elementwise<MAO_> = true,
                 enable_if_scalar<T_>        = true>
        _MATRIX_DIVIDED_BY_SCALAR_(MAO_ const&, T_);

        // Constructor (non-elementwise / scalar)
        template<typename MAO_                = MAO,
                 class    T_                  = T,
                 disable_if_elementwise<MAO_> = true,
                 enable_if_scalar<T_>         = true>
        _MATRIX_DIVIDED_BY_SCALAR_(MAO_ const&, T_);

        // Check if evaluation target is safe
        template<typename MAT>
        constexpr bool unsafe_evaluation_target(_MATRIX_<MAT> const&) const;

        // Check if expression references address of target matrix
        template<typename MAT>
        constexpr bool references_address(_MATRIX_<MAT> const*) const;

        // Evaluate the expression to the target matrix
        template<typename MAT>
        void eval(_MATRIX_<MAT>&) const;

        // Get the number of rows
        inline index_t rows() const;

        // Get the number of columns
        inline index_t cols() const;

        // Get the number of elements
        inline index_t size() const;

        // Get the element at the i-th row and j-th column
        dtype operator() (index_t, index_t) const;

        // Get the single-indexed element
        dtype operator() (index_t) const;
    };                                
}


#endif