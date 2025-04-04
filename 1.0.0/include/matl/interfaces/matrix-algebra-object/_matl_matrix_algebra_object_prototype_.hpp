#ifndef _MATL_MATRIX_ALGEBRA_OBJECT_PROTOTYPE__HPP
#define _MATL_MATRIX_ALGEBRA_OBJECT_PROTOTYPE__HPP

namespace matl
{
   /***********************************************************************************************
    * Matrix Algebra Object Prototype
    * 
    * This class is the base class for all matrix and matrix expression objects in the 'matl' 
    * library. It is used as an interface to specify which objects (aside from scalars) can be
    * used in matrix algebra operations.
    ***********************************************************************************************/
    template<typename MAO>
    class _MATRIX_ALGEBRA_OBJECT_
    {
      public:  

        // Matrix datatype
        typedef typename type_traits<MAO>::dtype dtype;

        // Row Getter
        index_t rows() const;

        // Coloumn Getter
        index_t cols() const;

        // Size Getter
        index_t size() const;
    };


}
#endif