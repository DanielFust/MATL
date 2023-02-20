//#ifndef _MATL_GENERIC_MAT_SCALAR_PROD__HPP
//#define _MATL_GENERIC_MAT_SCALAR_PROD__HPP
/**************************************************************************************************
 * Written by Daniel Fust
 **************************************************************************************************/
//#include "_matl_forward_declarations_.hpp"
//#include "_matl_type_traits_.hpp"
//#include "_matl_utils_.hpp"

namespace matl
{
   /***********************************************************************************************
    * class _GENERIC_MAT_SCALAR_PROD_
    * 
    * Template expression for the scalar multiplication of generic atomic matrix expressions. This
    * expression template is the default representation for matrix scalar multiplication. 
    * 
    * Dev Note: Non-atomic expression templates store row/column values as the objects they 
    * otherwise reference may be short lived (i.e. destroyed before this object).
    * 
    * Template Arguments:
    *     AE        (_ATOMIC_EXPR_ or Matrix) Atomic matrix expression or matrix
    *     T         Scalar data type
    ***********************************************************************************************/ 
    template<typename AE, class T>
    class _GENERIC_MAT_SCALAR_PROD_ : public _ATOMIC_EXPR_< _GENERIC_MAT_SCALAR_PROD_<AE,T> >
    {
      private:
        AE      const& mat_;
        T              val_;
        index_t        rows_;
        index_t        cols_;

        index_t block_size() const;

      public:
        // Return Type
        typedef typename promote<typename TypeTraits<AE>::dtype, T>::type dtype;

        // Constructor
        _GENERIC_MAT_SCALAR_PROD_(AE const& mat, T val);
        ~_GENERIC_MAT_SCALAR_PROD_(){std::cout<<"_GENERIC_MAT_SCALAR_PROD_ Destructor called"<<std::endl;}

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
//#endif