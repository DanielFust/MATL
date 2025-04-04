
namespace matl
{
   /**********************************************************************************************
    * class _GENERIC_MAT_TRANSPOSE_ : public _MATEXPRESSION_
    * 
    * Expression template for general matrix tranpositions. 
    * 
    * Template Parameters:
    *     MAT              Matrix or matrix expression 
    **********************************************************************************************/
    template<typename MAT>
    class _GENERIC_MAT_TRANSPOSE_ : public _MATEXPRESSION_< _GENERIC_MAT_TRANSPOSE_<MAT> >
    {
      private:
        MAT const& mat_;

        index_t block_size() const;
      
      public:

        // Return Type
        typedef typename TypeTraits<MAT>::dtype dtype;               

        // Constructor
        _GENERIC_MAT_TRANSPOSE_(MAT const& mat);
        //~_GENERIC_MAT_TRANSPOSE_(){}

        // Evaluate expression to target Matrix (generic)
        template<typename MAT_T>
        void eval(Matrix<MAT_T>& target) const;

        // Evaluate expression to new matrix (defaults to first matrix type)
        MAT eval() const;

        // Get indexed result of expression
        dtype operator() (index_t, index_t) const;

        // Number of rows in resultant matrix
        index_t rows() const;
        
        // Number of columns in resultant matrix
        index_t cols() const;

        // Total size of resultant matrix
        index_t size() const;

        // Check if there is memory conflict with target matrix
        template<typename MAT_T>
        bool memory_conflict(Matrix<MAT_T> const&) const;
    };

}
