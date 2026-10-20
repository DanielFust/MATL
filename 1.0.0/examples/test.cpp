#include <iostream>
#include "../include/matl.hpp"

template<typename T> class DummyExpression;
template<typename T> struct matl::type_traits<DummyExpression<T> >
{
    typedef T dtype;
    typedef matl::Matrix<T> matrix_type;
    typedef matl::Matrix<T> downcast_type;
};
template<typename T> struct matl::is_elementwise<DummyExpression<T> >
{
    static const bool value = false;
    typedef matl::Matrix<T> member_type;
};

template<typename T>
class DummyExpression : public matl::_EXPRESSION_<DummyExpression<T> >
{
  private:
    matl::index_t rows_;
    matl::index_t cols_;  
  public:
    typedef T dtype;
    DummyExpression(matl::index_t rows, matl::index_t cols) : rows_(rows), cols_(cols) {}

    matl::index_t rows() const {return rows_;}
    matl::index_t cols() const {return cols_;}
    matl::index_t size() const {return rows_*cols_;}
    template<typename MAT>
    void eval(matl::_MATRIX_<MAT>& mat) const
    {
        for (matl::index_t i = 1; i <= rows_; i++)
            for (matl::index_t j = 1; j <= cols_; j++)
                mat(i,j) = 1;
            
    }
    
};


template<typename T>
class MyClass
{
  public:

    //template<typename U = T, matl::enable_if_elementwise<U> = true> // doesn't work
    template<typename U = T, std::enable_if<matl::is_elementwise<U>::value>::type = true>
    MyClass(U& obj) {std::cout << "MyClass elementwise constructor" << std::endl;}


    template <std::enable_if<std::is_scalar<T>::value,bool>::type = true>
    MyClass(T& obj) {std::cout << "MyClass scalar constructor" << std::endl;}

    //MyClass(const MyClass<T>& other) {std::cout << "MyClass copy constructor" << std::endl;}

};

int main()
{
    using matl::RowMajorMatrix;
    using matl::FixedSizeRowMajorMatrix;
    using matl::Matrix;
    using matl::index_t;
    using matl::ScratchContainer;
    using matl::MemoryManager;
    typedef double dtype;
    
    index_t const m = 3;
    index_t const n = 3;
    dtype a = 1.0;
    dtype b = 2.0;

    matl::MemoryManager::ALLOCATE_MATRICES_IN_SCRATCH = false;

    /*
    Matrix<dtype> A(m,n);
    Matrix<dtype> B(m,n);
    Matrix<dtype> C(m,n);
    */
    
    
    Matrix<dtype,m,n> A;
    Matrix<dtype,m,n> B;
    Matrix<dtype,m,n> C;


    A = a;
    B = b;
    C = a;

    // Zero temporary matrices
    C = a + A*2 + B - b - B/2 + 4*C;

    // One temporary matrix
    C = A + B*A - B;
    
    // Zero temporary matrices
    C = A*B;

    // One temporary matrix
    C = A*C;
    
    // One temporary matrix
    C = B + A + transpose(A*B+C);

    // Zero temporary matrices
    C = transpose(A+B);
    
    // 
   
    std::cout << "Memory Manager Enabled: " << MemoryManager::get_memory_manager().enabled() << std::endl;
    std::cout << "Number of Scratch Containers: " << MemoryManager::get_memory_manager().containers().size() << std::endl;
    for (ScratchContainer& container : MemoryManager::get_memory_manager().containers())
    {
        std::cout << "Container size: " << container.size() << 
                     ", Bytes used: "   << container.bytes_used() 
                  << std::endl;
    }
    
    

    for (index_t i = 1; i <= m; i++)
    {
        for (index_t j = 1; j <= n; j++)
        {
            std::cout << C(i,j) << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}