#include <iostream>
#include <type_traits>

class A
{};
class B
{};
class C : public A
{};

template<typename T>
class Base
{};

class Derived : public Base<Derived>
{};


template<typename T> struct is_elementwise;

template<>
struct is_elementwise<A>
{
    static const bool value = true;
};

template<>
struct is_elementwise<B>
{
    static const bool value = false;
};

template<typename T>
struct is_elementwise<Base<T> >
{
    static const bool value = true;
};

void foo()
{
    A a;
    B b;
    C c;
    Derived d;

    auto a_val = is_elementwise<A>::value;
    auto b_val = is_elementwise<B>::value;

    std::cout << "Is A elementwise? " << a_val << std::endl;
    std::cout << "Is B elementwise? " << b_val << std::endl;

    std::enable_if<is_elementwise<A>::value, bool>::type c;
    std::enable_if<is_elementwise<B>::value, bool>::type d;
    std::enable_if<is_elementwise<C>::value, bool>::type c;
    std::enable_if<is_elementwise<Derived>::value, bool>::type d;
}

int main()
{
    //auto a = std::is_integral<int>::value;
    //auto b = std::is_integral<float>::value;

    //std::cout << "Is int integral? " << a << std::endl;
    //std::cout << "Is float integral? " << b << std::endl;

    //std::enable_if<std::is_scalar<T>::value, bool>::type
    //std::enable_if<false, float>::type c;

    foo();


}