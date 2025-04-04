#ifndef _MATL_ELEMENTWISE_PROTOTYPE__HPP
#define _MATL_ELEMENTWISE_PROTOTYPE__HPP
namespace matl
{
   /***********************************************************************************************
    * class _ELEMENTWISE_ 
    * 
    * Base class for elementwise operations on matrices.The class is used to define the interface
    * for elementwise operations and to provide a common interface for elementwise operations on
    * matrices and matrix expressions.
    * 
    * Template Parameters:
    *     EW         (_ELEMENTWISE_) Elementwise operation type
    **********************************************************************************************/
    template<typename EW>
    class _ELEMENTWISE_
    {
      public:

        // Resultant matrix data type (e.g. double, int, float, etc.)
        typedef typename type_traits<EW>::dtype dtype;

        // Get indexed result of expression
        dtype operator() (index_t, index_t) const;

        // Get single-indexed result of expression (called if same storage order)
        dtype operator() (index_t) const;

        // Get number of rows of resultant matrix
        index_t rows() const;

        // Get number of columns of resultant matrix
        index_t cols() const;

        // Get number of elements in resultant matrix
        index_t size() const;

    };
}
#endif