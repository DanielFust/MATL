
namespace matl
{
   /**********************************************************************************************
    * class _GENERIC_MAT_PROD_ : public _MATEXPRESSION_
    * 
    * Expression template for general matrix - matrix multiplications. 
    * 
    * Dev Note: Non-atomic expression templates store row/column values as the objects they 
    * otherwise reference may be short lived (i.e. destroyed before this object).
    * 
    * Template Parameters:
    *     MAT1              Matrix or atomic matrix expression
    *     MAT2              Matrix or atomic matrix expression     
    **********************************************************************************************/
    template<typename MAT1, typename MAT2>
    class _GENERIC_MAT_PROD_ : public _MATEXPRESSION_< _GENERIC_MAT_PROD_<MAT1,MAT2> >
    {
      private:
        MAT1    const& mat1_;
        MAT2    const& mat2_;
        index_t        rows_;
        index_t        cols_;

        index_t block_size() const;
      
      public:

        // Return Type
        typedef typename promote<typename TypeTraits<MAT1>::dtype, 
                                 typename TypeTraits<MAT2>::dtype >::type dtype;               

        // Constructor
        _GENERIC_MAT_PROD_(MAT1 const& mat1, MAT2 const& mat2);
        ~_GENERIC_MAT_PROD_(){std::cout<<"_GENERIC_MAT_PROD_ Destructor called"<<std::endl;}

        // Evaluate expression to target Matrix (generic)
        template<typename MAT_T>
        void eval(Matrix<MAT_T>& target) const;

        // Evaluate expression to new matrix (defaults to first matrix type)
        MAT1 eval() const;

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
