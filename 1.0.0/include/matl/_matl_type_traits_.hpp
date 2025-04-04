/**************************************************************************************************
 * File:   _matl_type_traits_.hpp
 * 
 * Description: Type traits for the MATL library.
 **************************************************************************************************/
#ifndef _MATL_TYPE_TRAITS__HPP
#define _MATL_TYPE_TRAITS__HPP

// Standard Library shortcuts/aliases
namespace matl
{
    // Enable if Scalar shorthand
    template<typename T>
    using enable_if_scalar = typename std::enable_if<std::is_scalar<T>::value, bool>::type;

    // Disable if Scalar shorthand
    template<typename T>
    using disable_if_scalar = typename std::enable_if<!std::is_scalar<T>::value, bool>::type;

    // Enable if Integral shorthand
    template<typename T>
    using enable_if_integral = typename std::enable_if<std::is_integral<T>::value, bool>::type;

    // Disable if Integral shorthand
    template<typename T>
    using disable_if_integral = typename std::enable_if<!std::is_integral<T>::value, bool>::type;


    /******************* Type Promotion  *********************/

    // Type promotion struct
    template<typename... TYPES> struct promote;

    // End of type promotion recursion
    template<>  
    struct promote<>
    {typedef bool type;};

    // General purpose type promotion
    template<typename FIRST_TYPE, typename... OTHER_TYPES>
    struct promote<FIRST_TYPE, OTHER_TYPES...> 
    {typedef decltype(FIRST_TYPE() + typename promote<OTHER_TYPES...>::type()) type;};

}

// MATL Type Traits declarations
#include "./type-traits/_matl_type_traits_declarations_.hpp"

namespace matl
{


    /************ Elementwise Enable Disable shortcuts *************/

    // Enable if elementwise shorthand
    template<typename T>
    using enable_if_elementwise = typename std::enable_if<is_elementwise<T>::value, bool>::type;

    // Disable if elementwise shorthand
    template<typename T>
    using disable_if_elementwise = typename std::enable_if<!is_elementwise<T>::value, bool>::type;



    /******************* Default Matrix Types  *********************/

    // Default matrix type selector
    template<class T, index_t ... DIMS> struct default_matrix_selector;

    // Default matrix type selector specialization
    template<class T>
    struct default_matrix_selector<T>
    {
        typedef RowMajorMatrix<T> type;
    };

    // Default matrix type selector specialization
    template<class T, index_t M, index_t N>
    struct default_matrix_selector<T, M, N>
    {
        typedef FixedSizeRowMajorMatrix<T, M, N> type;
    };


    // Alias for default dynamic and fixed-size matrix types
    template<class T, index_t ... DIMS>
    using Matrix = typename default_matrix_selector<T,DIMS...>::type;



    /****************** unqualified type casting ******************/
    // Casts down to the unqualified type name of the MATL object even for multiple
    // levels of inheritance.

    // non-const specialization
    template<typename T>
    typename type_traits<T>::unqualified_type& unqualified_cast(T& obj)
    {return static_cast<type_traits<T>::unqualified_type&>(obj);}

    // const specialization
    template<typename T>
    typename type_traits<T>::unqualified_type const& unqualified_cast(T const& obj)
    {return static_cast<type_traits<T>::unqualified_type const&>(obj);}



    /*************** Concepts *******************/

    // Scalar
    template<class T>
    concept Scalar = std::is_scalar<T>::value;

    // Integer
    template<class T>
    concept Integer = std::is_integral<T>::value;

    // FixedSize
    template<typename T>
    concept FixedSize = is_fixed_size<T>::value;

    // MatrixAlgebraObject
    template<typename T>
    concept MatrixAlgebraObject = is_matrix_algebra_object<T>::value;

    // Elementwise
    template<typename T>
    concept Elementwise = is_elementwise<T>::value && is_matrix_algebra_object<T>::value;

    // NonElementwise
    template<typename T>
    concept NonElementwise = !is_elementwise<T>::value && is_matrix_algebra_object<T>::value;

} 


/******************* Specializations  *********************/
#include "./type-traits/_matl_type_traits_specializations_.hpp"


#endif