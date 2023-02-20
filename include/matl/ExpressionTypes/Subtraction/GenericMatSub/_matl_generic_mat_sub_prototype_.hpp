//#ifndef _MATL_GENERIC_MAT_SUB__HPP
//#define _MATL_GENERIC_MAT_SUB__HPP
/**************************************************************************************************
 * Written by Daniel Fust
 **************************************************************************************************/
//#include "_matl_forward_declarations_.hpp"
//#include "_matl_type_traits_.hpp"
//#include "_matl_utils_.hpp"

namespace matl
{
   /***********************************************************************************************
    * class _GENERIC_MAT_SUB_
    * 
    * Template expression for the subtraction of generic atomic matrix expressions. This expression
    * template is the default representation for matrix subtraction. This expression template 
    * supports lazy execution and some cache blocking, however more specific types are recommended
    * to take advantage of specific data storage architectures of Matrix types.
    * 
    * Template Arguments:
    *     AE_T1     (_ATOMIC_EXPR_) Atomic matrix expression
    *     AE_T2     (_ATOMIC_EXPR_) Atomic matrix expression
    ***********************************************************************************************/ 
    template<typename AE_T1, typename AE_T2>
    class _GENERIC_MAT_SUB_ : public _ATOMIC_EXPR_< _GENERIC_MAT_SUB_<AE_T1,AE_T2> >
    {
      private:
        AE_T1 const& mat1_;
        AE_T2 const& mat2_;

        index_t block_size() const;

      public:
        // Return Type
        typedef typename promote<typename TypeTraits<AE_T1>::dtype, 
                                 typename TypeTraits<AE_T2>::dtype>::type dtype;

        // Constructor
        _GENERIC_MAT_SUB_(AE_T1 const& mat1, AE_T2 const& mat2);

        // Compute and return (double) indexed value
        dtype operator() (index_t, index_t) const;

        // Evaluate expression 
        //Matrix<typename TypeTraits<AE_T1>::mat_type> eval() const;

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
//#endif