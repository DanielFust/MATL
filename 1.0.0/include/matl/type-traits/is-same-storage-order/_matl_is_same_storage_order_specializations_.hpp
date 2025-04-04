#ifndef _MATL_IS_SAME_STORAGE_ORDER__HPP
#define _MATL_IS_SAME_STORAGE_ORDER__HPP

namespace matl 
{
    /***********************************************************************************************
     * is_same_storage_order
     * 
     * Type trait to check if two matrix or expression types have the same storage order (row-major 
     * or column-major).
     **********************************************************************************************/
    
    // Recusion
    template<typename FIRST_TYPE, typename ... OTHER_TYPES>
    struct is_same_storage_order
    {
        static constexpr bool value = (is_rowmajor_storage<FIRST_TYPE>::value && is_rowmajor_storage<OTHER_TYPES...>::value) ||
                                      (is_colmajor_storage<FIRST_TYPE>::value && is_colmajor_storage<OTHER_TYPES...>::value);                            
    };
    
    // End of recusion (single type just needs to have a storage order)
    template<typename T>
    struct is_same_storage_order<T>
    {static constexpr bool value = is_rowmajor_storage<T>::value || is_colmajor_storage<T>::value;};
}
#endif