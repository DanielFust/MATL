/**************************************************************************************************
 * File: _matl_owned_implementation_.hpp
 *
 * Description: Implementation of the owned matrix expression class.
 **************************************************************************************************/
#ifndef MATL_OWNED_IMPLEMENTATION_HPP
#define MATL_OWNED_IMPLEMENTATION_HPP
namespace matl::expr
{
   /***********************************************************************************************
    * Constructor
    *
    * Takes ownership of a matrix. The matrix's data is moved, not copied.
    *
    * Parameters:
    *     mat          Matrix to own
    ***********************************************************************************************/
    template<typename MAT>
    Owned<MAT>::Owned(MAT mat) :
        mat_(std::make_shared<const MAT>(std::move(mat)))
    {}

   /***********************************************************************************************
    * MAT const& matrix() const
    *
    * Returns:
    *     The owned matrix
    ***********************************************************************************************/
    template<typename MAT>
    MAT const& Owned<MAT>::matrix() const
    {return *this->mat_;}

   /***********************************************************************************************
    * index_t rows() const
    *
    * Returns:
    *      index_t      Number of rows in the owned matrix
    **********************************************************************************************/
    template<typename MAT>
    index_t Owned<MAT>::rows() const
    {return this->mat_->rows();}

   /***********************************************************************************************
    * index_t cols() const
    *
    * Returns:
    *     index_t      Number of columns in the owned matrix
    ***********************************************************************************************/
    template<typename MAT>
    index_t Owned<MAT>::cols() const
    {return this->mat_->cols();}

   /***********************************************************************************************
    * index_t size() const
    *
    * Returns:
    *     index_t      Number of elements in the owned matrix
    ***********************************************************************************************/
    template<typename MAT>
    index_t Owned<MAT>::size() const
    {return this->mat_->size();}

   /***********************************************************************************************
    * bool unsafe_evaluation_target(TARGET const&)
    *
    * The owned matrix is a separate object from any evaluation target, so evaluation is always
    * safe.
    *
    * Returns:
    *     'false'
    ***********************************************************************************************/
    template<typename MAT>
    template<MatrixType TARGET>
    bool Owned<MAT>::unsafe_evaluation_target(TARGET const&) const
    {return false;}

   /***********************************************************************************************
    * bool references_address(void const*)
    *
    * The owned matrix is a separate object from any evaluation target.
    *
    * Returns:
    *     'false'
    ***********************************************************************************************/
    template<typename MAT>
    bool Owned<MAT>::references_address(void const*) const
    {return false;}

   /***********************************************************************************************
    * void eval(TARGET& target) const
    *
    * Copy the owned matrix to the target matrix.
    *
    * Parameters:
    *     target       Target matrix
    ***********************************************************************************************/
    template<typename MAT>
    template<MatrixType TARGET>
    void Owned<MAT>::eval(TARGET& target) const
    {detail::assign_elementwise(target, *this);}
}

namespace matl::detail
{
   /***********************************************************************************************
    * struct evaluator<expr::Owned>
    *
    * Reads elements of the owned matrix directly.
    ***********************************************************************************************/
    template<typename MAT>
    struct evaluator< expr::Owned<MAT> >
    {
        typedef typename expr::Owned<MAT>::dtype dtype;

        evaluator<MAT> mat_;

        explicit evaluator(expr::Owned<MAT> const& expr) : mat_(expr.matrix()) {}

        dtype operator() (index_t i, index_t j) const {return mat_(i,j);}
        dtype operator() (index_t i)            const {return mat_(i);}
    };
}
#endif
