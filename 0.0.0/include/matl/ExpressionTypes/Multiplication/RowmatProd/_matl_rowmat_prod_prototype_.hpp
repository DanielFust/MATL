//#ifndef _MATL_ROWMAT_PROD__HPP
//#define _MATL_ROWMAT_PROD__HPP

//#include "_matl_forward_declarations_.hpp"
//#include "_matl_type_traits_.hpp"


namespace matl
{
   /**********************************************************************************************
    * class _ROWMAT_PROD_ : public _MATEXPRESSION_
    * 
    * Expression template for matrix products between two dense row-major matricies (RowMatrix).
    * 
    * Template Parameters:
    *     T1 (implicit)        Data type of first RowMatrix
    *     T2 (implicit)        Data type of second RowMatrix     
    **********************************************************************************************/
    template<typename RMAE1, typename RMAE2>
    class _ROWMAT_PROD_ : public _MATEXPRESSION_< _ROWMAT_PROD_<RMAE1,RMAE2> >
    {
      private:
        RowMatrix<typename promote<typename TypeTraits<RMAE1>::dtype, 
                                   typename TypeTraits<RMAE2>::dtype >::type> mat_;
      public:

        // Return Type
        typedef typename promote<typename TypeTraits<RMAE1>::dtype, 
                                 typename TypeTraits<RMAE2>::dtype >::type dtype;

        // Constructor
        _ROWMAT_PROD_(RMAE1 const& mat1, RMAE2 const& mat2);

        // Evaluate expression 
        RowMatrix<dtype> eval() const;

        // Evaluate expression to target Matrix (generic)
        template<typename MAT_T>
        void eval(Matrix<MAT_T>& target) const;

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

//#endif