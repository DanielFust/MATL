#ifndef _MATL_FIXED_SIZE_PROTOTYPE__HPP
#define _MATL_FIXED_SIZE_PROTOTYPE__HPP

namespace matl 
{
   /***********************************************************************************************
    * class _FIXED_SIZE_
    * 
    * Interface class for fixed-size matrices and expressions. All fixed-size Matrix Algebra 
    * objects must be able to provide their dimensions and size at compile time.
    **********************************************************************************************/ 
    template<typename FS>
    class _FIXED_SIZE_ : public _MATRIX_ALGEBRA_OBJECT_< _FIXED_SIZE_<FS> >
    {
      public:  

        // Data type
        typedef typename type_traits<FS>::dtype dtype;

        // Row Getter
        static constexpr index_t ROWS = FS::ROWS;

        // Coloumn Getter
        static constexpr index_t COLS = FS::COLS;

        // Size Getter
        static constexpr index_t SIZE = FS::SIZE;
    };
}

#endif