#include <iostream>
#include <type_traits>
void modify_int(int const& number)
{
    //int* num_ptr = static_cast<int*>(&number);      // static_cast can't cast away const
    //int* num_ptr = reinterpret_cast<int*>(&number); // reinterpret_cast can't cast away const

    //int* num_ptr = (int*) &number; // Works!
    int* num_ptr = const_cast<int*>(&number); // Also works!
}

template <class T, std::enable_if<std::is_integral<T>::value> >
void foo(T& value)
{
    std::cout << "Integral value: " << value << std::endl;
}

int main()
{
    auto a = std::is_integral<int>::value;
    auto b = std::is_integral<float>::value;

    typedef int const& constant_integer;

    constant_integer a_int = 5;
    std::cout << "Before const_cast: " << a_int << std::endl;
    const_cast<int&>(a_int) += 5;
    std::cout << "After const_cast: " << a_int << std::endl;

    //std::cout << "Is int integral? " << a << std::endl;
    //std::cout << "Is float integral? " << b << std::endl;

}


