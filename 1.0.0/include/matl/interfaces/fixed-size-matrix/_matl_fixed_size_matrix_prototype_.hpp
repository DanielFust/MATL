#ifndef _MATL_FIXED_SIZE_MATRIX_PROTOTYPE__HPP
#define _MATL_FIXED_SIZE_MATRIX_PROTOTYPE__HPP
namespace matl
{
   /***********************************************************************************************
    * class _FIXED_SIZE_MATRIX_
    * 
    * This class is the base class for all fixed-size matrices in the 'matl' library. Matrix
    * types that inherit from this interface class may be used in fixed-size expressions that
    * do not require dynamic memory allocation.
    ***********************************************************************************************/
    template<typename FSM>
    class _FIXED_SIZE_MATRIX_ : public _MATRIX_< _FIXED_SIZE_MATRIX_<FSM> >,
                                public _FIXED_SIZE_< _FIXED_SIZE_MATRIX_<FSM> >
    {
      public:  

        /******************* _MATRIX_ interface *******************/
        // Matrix datatype
        typedef typename type_traits<FSM>::dtype dtype;

        // Matrix copy
        template<typename MAT>
        FSM& operator= (_MATRIX_<MAT>& mat);

        // Expression evaluation trigger
        template<typename EXPR>
        FSM& operator=(_EXPRESSION_<EXPR> const& expr);

        // Single index getter
        dtype operator() (index_t i) const;

        // Double index getter
        dtype operator() (index_t i, index_t j) const;

        // Single index accessor
        dtype& operator() (index_t i);

        // Double index accessor
        dtype& operator() (index_t i, index_t j);

        // Row Getter
        index_t rows() const;

        // Coloumn Getter
        index_t cols() const;

        // Size Getter
        index_t size() const;

        /***************** _FIXED_SIZE_ interface ******************/

        // Static row getter
        static constexpr index_t ROWS = FSM::ROWS;

        // Static column getter
        static constexpr index_t COLS = FSM::COLS;

        // Static size getter
        static constexpr index_t SIZE = FSM::SIZE;
    };
}
#endif    