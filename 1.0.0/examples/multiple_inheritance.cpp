#include <iostream>

template<typename T>
class BaseA
{
  public:  
    void hello()
    {static_cast<T*>(this)->hello();} 

    void diamond()
    {static_cast<T*>(this)->diamond();}
};

template<typename T>
class BaseB
{
  public:
    void goodbye()
    {static_cast<T*>(this)->goodbye();}
    
    void diamond()
    {static_cast<T*>(this)->diamond();}
};


class DerivedA : public BaseA<DerivedA>,
                 public BaseB<DerivedA>
{
  public:  
    void hello()
    {std::cout << "Hello from DerivedA!" << std::endl;}

    void goodbye()
    {std::cout << "Goodbye from DerivedA!" << std::endl;}

    void diamond()
    {std::cout << "Diamond from DerivedA!" << std::endl;}
};


template<typename T>
void hello(BaseA<T>& b)
{
    b.hello();
    b.diamond();
}

template<typename T>
void goodbye(BaseB<T>& b)
{
    b.goodbye();
    b.diamond();
}


int main()
{
    DerivedA a;
    a.hello();
    hello(a);

    a.goodbye();
    goodbye(a);
    return 0;
}