/**************************************************************************************************
 * File: _matl_matrix_times_matrix_implementation_.hpp
 *
 * Description: Implementations for the matrix times matrix expression template.
 **************************************************************************************************/
#ifndef _MATL_MATRIX_TIMES_MATRIX_IMPLEMENTATION__HPP
#define _MATL_MATRIX_TIMES_MATRIX_IMPLEMENTATION__HPP
namespace matl
{
   /***********************************************************************************************
    * Constructor
    *
    * Stores the operands: expressions are copied, matrices are referenced. Nothing is computed.
    *
    * Parameters:
    *     mat1         Left operand (matrix or expression)
    *     mat2         Right operand (matrix or expression)
    ***********************************************************************************************/
    template<typename MAO1, typename MAO2>
    template<typename OP1, typename OP2>
    _MATRIX_TIMES_MATRIX_<MAO1,MAO2>::_MATRIX_TIMES_MATRIX_(OP1 const& mat1, OP2 const& mat2) :
        mat1_(unqualified_cast(mat1)),
        mat2_(unqualified_cast(mat2))
    {}

   /***********************************************************************************************
    * lhs_type const& lhs() const
    *
    * Returns:
    *     Left operand
    ***********************************************************************************************/
    template<typename MAO1, typename MAO2>
    typename _MATRIX_TIMES_MATRIX_<MAO1,MAO2>::lhs_type const& _MATRIX_TIMES_MATRIX_<MAO1,MAO2>::lhs() const
    {return this->mat1_;}

   /***********************************************************************************************
    * rhs_type const& rhs() const
    *
    * Returns:
    *     Right operand
    ***********************************************************************************************/
    template<typename MAO1, typename MAO2>
    typename _MATRIX_TIMES_MATRIX_<MAO1,MAO2>::rhs_type const& _MATRIX_TIMES_MATRIX_<MAO1,MAO2>::rhs() const
    {return this->mat2_;}

   /***********************************************************************************************
    * index_t rows() const
    *
    * Returns:
    *     index_t      Number of rows in resultant matrix
    ***********************************************************************************************/
    template<typename MAO1, typename MAO2>
    index_t _MATRIX_TIMES_MATRIX_<MAO1,MAO2>::rows() const
    {return mat1_.rows();}

   /***********************************************************************************************
    * index_t cols() const
    *
    * Returns:
    *     index_t      Number of columns in resultant matrix
    ***********************************************************************************************/
    template<typename MAO1, typename MAO2>
    index_t _MATRIX_TIMES_MATRIX_<MAO1,MAO2>::cols() const
    {return mat2_.cols();}

   /***********************************************************************************************
    * index_t size() const
    *
    * Returns:
    *     index_t      Number of elements in resultant matrix
    ***********************************************************************************************/
    template<typename MAO1, typename MAO2>
    index_t _MATRIX_TIMES_MATRIX_<MAO1,MAO2>::size() const
    {return mat1_.rows() * mat2_.cols();}

   /***********************************************************************************************
    * bool unsafe_evaluation_target(_MATRIX_ const&)
    *
    * Check if computing the product directly into the target matrix would corrupt memory.
    *
    * The product reads matrix operands directly while it writes the target, so it is unsafe if
    * either operand IS the target (e.g. C = A*C). Expression operands are evaluated into
    * temporaries before any element of the target is written, so they are always safe, even if
    * they read the target (e.g. C = (A*C)*B or C = transpose(C)*B).
    *
    * Parameters:
    *     target       Target matrix
    * Returns:
    *     'true' if computing the product directly into the target would corrupt memory
    *     'false' otherwise
    ***********************************************************************************************/
    template<typename MAO1, typename MAO2>
    template<typename MAT>
    constexpr bool _MATRIX_TIMES_MATRIX_<MAO1,MAO2>::unsafe_evaluation_target(_MATRIX_<MAT> const& target) const
    {
        bool unsafe = false;
        if constexpr (MatrixType<lhs_type>) unsafe = unsafe || matl::references_address(&target, this->mat1_);
        if constexpr (MatrixType<rhs_type>) unsafe = unsafe || matl::references_address(&target, this->mat2_);
        return unsafe;
    }

   /***********************************************************************************************
    * bool references_address(_MATRIX_ const*)
    *
    * Check if the expression references the target matrix.
    *
    * Parameters:
    *     target       Target matrix
    * Returns:
    *     'true' if the expression references the target matrix, 'false' otherwise
    ***********************************************************************************************/
    template<typename MAO1, typename MAO2>
    template<typename MAT>
    constexpr bool _MATRIX_TIMES_MATRIX_<MAO1,MAO2>::references_address(_MATRIX_<MAT> const* target) const
    {return (matl::references_address(target,this->mat1_) || matl::references_address(target,this->mat2_));}

   /***********************************************************************************************
    * void compute_product(_MATRIX_<MAT>& target, _MATRIX_TIMES_MATRIX_ const& product)
    *
    * Computes a matrix product into the target matrix. Matrix operands are read directly;
    * expression operands are first evaluated into temporary matrices, which are released as soon
    * as the product has been computed.
    *
    * Parameters:
    *     target       Target matrix (must not be read by the product, see unsafe_evaluation_target)
    *     product      Product expression
    ***********************************************************************************************/
    template<typename MAT, typename MAO1, typename MAO2>
    void compute_product(_MATRIX_<MAT>& target, _MATRIX_TIMES_MATRIX_<MAO1,MAO2> const& product)
    {
        typedef _MATRIX_TIMES_MATRIX_<MAO1,MAO2> expression_type;
        typedef typename expression_type::dtype  dtype;

        const product_operand_t<typename expression_type::lhs_type> lhs(product.lhs());
        const product_operand_t<typename expression_type::rhs_type> rhs(product.rhs());

        const index_t M = product.rows();
        const index_t N = product.cols();
        const index_t K = product.lhs().cols();

        for (index_t i = 1; i <= M; i++)
        {
            for (index_t j = 1; j <= N; j++)
            {
                dtype sum = 0;
                for (index_t k = 1; k <= K; k++)
                    sum += lhs(i,k) * rhs(k,j);
                target(i,j) = sum;
            }
        }
    }

   /***********************************************************************************************
    * void eval(_MATRIX_<MAT>& target) const
    *
    * Computes the product directly into the target matrix.
    *
    * Parameters:
    *     target       Target matrix (must be a safe evaluation target with matching dimensions)
    ***********************************************************************************************/
    template<typename MAO1, typename MAO2>
    template<typename MAT>
    void _MATRIX_TIMES_MATRIX_<MAO1,MAO2>::eval(_MATRIX_<MAT>& target) const
    {compute_product(target, *this);}

   /***********************************************************************************************
    * struct evaluator<_MATRIX_TIMES_MATRIX_>
    *
    * Evaluator for a product nested inside a larger expression. The product is computed into a
    * temporary matrix when the evaluator is constructed (before any element of the assignment
    * target is written), and the temporary is released when the evaluator is destroyed at the
    * end of the assignment.
    ***********************************************************************************************/
    template<typename MAO1, typename MAO2>
    struct evaluator< _MATRIX_TIMES_MATRIX_<MAO1,MAO2> >
    {
        typedef _MATRIX_TIMES_MATRIX_<MAO1,MAO2>                    expression_type;
        typedef typename expression_type::dtype                     dtype;
        typedef typename is_fixed_size<expression_type>::matrix_type result_type;

        result_type result_;

        explicit evaluator(expression_type const& expr) :
            result_(temporary_matrix<result_type>(expr.rows(), expr.cols()))
        {compute_product(result_, expr);}

        dtype operator() (index_t i, index_t j) const {return result_(i,j);}
        dtype operator() (index_t i)            const {return result_(i);}
    };
}
#endif
