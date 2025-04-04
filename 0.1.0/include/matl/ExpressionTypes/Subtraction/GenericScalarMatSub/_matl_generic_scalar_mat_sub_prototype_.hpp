
namespace matl
{
   /***********************************************************************************************
    * class _GENERIC_MAT_SCALAR_SUB_
    * 
    * Template expression for the subtraction of a scalar from a generic matrix or atomic matrix 
    * expression.
    * 
    * Template Arguments:
    *     AE_T1     (_ATOMIC_EXPR_) Atomic matrix expression
    *     AE_T2     (_ATOMIC_EXPR_) Atomic matrix expression
    ***********************************************************************************************/ 
    template<class T, typename EXPR>
    class _GENERIC_SCALAR_MAT_SUB_ : public _MATEXPRESSION_< _GENERIC_SCALAR_MAT_SUB_<T,EXPR> >
    {
      private:
      
        EXPR const& mat_;
        T           val_;

        index_t block_size() const;

      public:

        // Return Type
        typedef typename TypeTraits<EXPR>::dtype dtype;

        // Constructor
        _GENERIC_SCALAR_MAT_SUB_(T val, EXPR const& mat);

        // Compute and return (double) indexed value
        dtype operator() (index_t, index_t) const;

        // Evaluate expression to target Matrix (generic)
        template<typename MAT_T>
        void eval(Matrix<MAT_T>& target) const;

        // Number of rows in resultant matrix
        index_t rows() const;
        
        // Number of columns in resultant matrix
        index_t cols() const;

        // Total size of resultant matrix
        index_t size() const;

    };

}
