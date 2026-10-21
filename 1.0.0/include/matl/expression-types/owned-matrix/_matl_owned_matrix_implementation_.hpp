/**************************************************************************************************
 * File: _matl_owned_matrix_implementation_.hpp
 *
 * Description: Implementations for the owned matrix expression class.
 **************************************************************************************************/
#ifndef _MATL_OWNED_MATRIX_IMPLEMENTATION__HPP
#define _MATL_OWNED_MATRIX_IMPLEMENTATION__HPP
namespace matl
{
   /***********************************************************************************************
    * Constructor
    *
    * Takes ownership of a matrix. The matrix's data is moved, not copied (for matrices borrowing
    * scratch memory, the ScratchContainer pairing moves with it).
    *
    * Parameters:
    *     mat          Matrix to own
    ***********************************************************************************************/
    template<typename MAT>
    _OWNED_MATRIX_<MAT>::_OWNED_MATRIX_(MAT mat) :
        mat_(std::make_shared<const MAT>(std::move(mat)))
    {}

   /***********************************************************************************************
    * MAT const& matrix() const
    * 
    * Returns:
    *     The owned matrix
    ***********************************************************************************************/
    template<typename MAT>
    MAT const& _OWNED_MATRIX_<MAT>::matrix() const
    {return *this->mat_;}

   /***********************************************************************************************
    * index_t rows() const
    *
    * Returns:
    *      index_t      Number of rows in the owned matrix
    **********************************************************************************************/
    template<typename MAT>
    index_t _OWNED_MATRIX_<MAT>::rows() const
    {return this->mat_->rows();}

   /***********************************************************************************************
    * index_t cols() const
    *
    * Returns:
    *     index_t      Number of columns in the owned matrix
    ***********************************************************************************************/
    template<typename MAT>
    index_t _OWNED_MATRIX_<MAT>::cols() const
    {return this->mat_->cols();}

   /***********************************************************************************************
    * index_t size() const
    *
    * Returns:
    *     index_t      Number of elements in the owned matrix
    ***********************************************************************************************/
    template<typename MAT>
    index_t _OWNED_MATRIX_<MAT>::size() const
    {return this->mat_->size();}

   /***********************************************************************************************
    * bool unsafe_evaluation_target(_MATRIX_ const&)
    *
    * The owned matrix is a separate object from any evaluation target, so evaluation is always
    * safe.
    *
    * Returns:
    *     'false'
    ***********************************************************************************************/
    template<typename MAT>
    template<typename MAT2>
    constexpr bool _OWNED_MATRIX_<MAT>::unsafe_evaluation_target(_MATRIX_<MAT2> const&) const
    {return false;}

   /***********************************************************************************************
    * bool references_address(_MATRIX_ const*)
    *
    * The owned matrix is a separate object from any evaluation target.
    *
    * Returns:
    *     'false'
    ***********************************************************************************************/
    template<typename MAT>
    template<typename MAT2>
    constexpr bool _OWNED_MATRIX_<MAT>::references_address(_MATRIX_<MAT2> const*) const
    {return false;}

   /***********************************************************************************************
    * void eval(_MATRIX_<MAT2>& target) const
    *
    * Copy the owned matrix to the target matrix.
    *
    * Parameters:
    *     target       Target matrix
    ***********************************************************************************************/
    template<typename MAT>
    template<typename MAT2>
    void _OWNED_MATRIX_<MAT>::eval(_MATRIX_<MAT2>& target) const
    {assign_elementwise(target, *this);}

   /***********************************************************************************************
    * struct evaluator<_OWNED_MATRIX_>
    * 
    * Reads elements of the owned matrix directly.
    ***********************************************************************************************/
    template<typename MAT>
    struct evaluator< _OWNED_MATRIX_<MAT> >
    {
        typedef typename _OWNED_MATRIX_<MAT>::dtype dtype;

        evaluator<MAT> mat_;

        explicit evaluator(_OWNED_MATRIX_<MAT> const& expr) : mat_(expr.matrix()) {}

        dtype operator() (index_t i, index_t j) const {return mat_(i,j);}
        dtype operator() (index_t i)            const {return mat_(i);}
    };
}
#endif
