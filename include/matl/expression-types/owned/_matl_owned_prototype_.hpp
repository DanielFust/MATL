/**************************************************************************************************
 * File: _matl_owned_prototype_.hpp
 *
 * Description: Prototype for the owned matrix expression class.
 **************************************************************************************************/
#ifndef MATL_OWNED_PROTOTYPE_HPP
#define MATL_OWNED_PROTOTYPE_HPP

namespace matl
{
    namespace detail
    {
        // An owned matrix has the traits of the matrix it owns
        template<typename MAT>
        struct traits< expr::Owned<MAT> > : traits<MAT> {};

        template<typename MAT>
        inline constexpr bool is_owned_v< expr::Owned<MAT> > = true;
    }

    namespace expr
    {
       /*******************************************************************************************
        * class Owned<MAT>
        *
        * Elementwise expression that owns a matrix. Temporary (rvalue) matrices used as operands,
        * e.g. 'A + eye<double>(3)', are moved into an Owned expression so that the expression
        * remains valid after the temporary would otherwise have been destroyed.
        *
        * The matrix is held through a shared pointer to const, so copies of the expression (and
        * of expressions containing it) share the matrix instead of copying its data. The matrix
        * is never modified after it is moved in.
        *
        * Template Parameters:
        *     MAT        Matrix type (e.g. RowMajorMatrix<T>, FixedSizeRowMajorMatrix<T,M,N>)
        *******************************************************************************************/
        template<typename MAT>
        class Owned : public ExpressionBase< Owned<MAT> >
        {
          public:
            // Resultant matrix data type (e.g. double, int, float, etc.)
            typedef detail::dtype_t<MAT> dtype;

          private:
            // Owned matrix
            std::shared_ptr<const MAT> mat_;

          public:

            // Constructor (takes ownership of the matrix)
            explicit Owned(MAT mat);

            // The owned matrix
            MAT const& matrix() const;

            // Check if evaluation target is safe
            template<WritableMatrix TARGET>
            bool unsafe_evaluation_target(TARGET const&) const;

            // Check if the expression reads memory that may hold elements of a target (see detail::Footprint)
            bool reads(detail::Footprint const&) const;

            // Evaluate the expression to the target matrix
            template<WritableMatrix TARGET>
            void eval(TARGET&) const;

            // Get the number of rows
            index_t rows() const;

            // Get the number of columns
            index_t cols() const;

            // Get the number of elements
            index_t size() const;
        };
    }
}

#endif
