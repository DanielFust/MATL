#include <iostream>
//#include "./include/matl.hpp"




class MyClass1
{
    void my_method()
    {
        std::cout << "Hi, I'm an instance of MyClass1" << std::endl;
    }
};

class MyClass2
{
    void my_method()
    {
        std::cout << "Hi, I'm an instance of MyClass2" << std::endl;
    }
};



class Node
{
  private:  
    void* data_;
    Node* next_;
  public:
    template<class T>  
    Node(T* data) : data_((void*)data), next_(nullptr) {}

    Node(void* data) : data_(data), next_(nullptr) {}

    void* data()
    {return this->data_;}
};


class ExpressionManager
{
  private:
    Node* head_;
  public:  
    ExpressionManager() : head_(nullptr){}
    
};

template<class T>
struct Container
{
    T* val_;
};




 


int main()
{
    MyClass1 obj1;
    MyClass2 obj2;

    Node node1((void*)&obj1);
    Node node2(&obj2);
    
    auto obj3 = node1.data();
    auto obj4 = node1.data();
    

    return 0;
}


