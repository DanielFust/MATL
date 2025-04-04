#include <iostream>

template<typename CLASS_T>
class BaseClass {};

template<class DATA_T>
class DerivedType1 : public BaseClass< DerivedType1<DATA_T> > {};

template<class DATA_T>
class DerivedType2 : public BaseClass< DerivedType2<DATA_T> > {};

// Operator template for DerivedType - DerivedType multiplications
template<typename CLASS_T1, typename CLASS_T2>
int operator* (BaseClass<CLASS_T1> const& obj1, BaseClass<CLASS_T2> const& obj2)
{
    // Multiplication operation
    return 0; // would actually return a different class object, though this isn't relevant to the issue
}


/*
// C++20 only
template<typename T>
concept Scalar = std::is_integral_v<T> || std::is_floating_point_v<T>;

template<typename CLASS_T, Scalar T>
int operator* (BaseClass<CLASS_T> const& obj, T scalar)
{
    // Multiplication operation
    return 1;
}
*/


template<typename T>
using enable_if_scalar = std::enable_if_t<std::is_scalar_v<T>, bool>;

template<typename CLASS_T, typename T, enable_if_scalar<T> = true>
int operator* (BaseClass<CLASS_T> const& obj, T scalar) 
{
    // Multiplication operation
    return 1;
}


// Main script
int main()
{
    DerivedType1<float> obj1;
    DerivedType2<float> obj2;

    int val1 = obj1*obj2;
    int val2 = obj2*2.0;

    std::cout << val1 << std::endl;
    std::cout << val2 << std::endl;

    return 0;
}
