#include <iostream>
#include "../include/matl.hpp"

// A user-defined expression type. It describes itself to MATL with a detail::traits
// specialization and derives from ExpressionBase, which makes it an ExpressionType that can be
// assigned to matrices and combined with other operands.
template<typename T> class DummyExpression;
template<typename T> struct matl::detail::traits<DummyExpression<T> >
{
    typedef T dtype;
    static constexpr bool          fixed_size    = false;
    static constexpr matl::index_t rows          = 0;
    static constexpr matl::index_t cols          = 0;
    static constexpr matl::detail::StorageOrder storage_order = matl::detail::StorageOrder::RowMajor;
};

template<typename T>
class DummyExpression : public matl::ExpressionBase<DummyExpression<T> >
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

    // Reads no matrices, so it is safe to evaluate into any target
    template<matl::MatrixType MAT>
    bool unsafe_evaluation_target(MAT const&) const {return false;}
    bool references_address(void const*) const {return false;}

    template<matl::MatrixType MAT>
    void eval(MAT& mat) const
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

    // Constructor for matrices and expressions
    template<matl::MatrixAlgebraObject U>
    MyClass(U& /*obj*/) {std::cout << "MyClass matrix/expression constructor" << std::endl;}

    // Constructor for scalars
    MyClass(T& /*obj*/) requires matl::Scalar<T> {std::cout << "MyClass scalar constructor" << std::endl;}

    //MyClass(const MyClass<T>& other) {std::cout << "MyClass copy constructor" << std::endl;}

};

int main()
{
    using matl::RowMajorMatrix;
    using matl::FixedSizeRowMajorMatrix;
    using matl::Matrix;
    using matl::index_t;
    using matl::MemoryManager;
    typedef double dtype;
    
    index_t const m = 3;
    index_t const n = 3;
    dtype a = 1.0;
    dtype b = 2.0;

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

    // User-defined expression type
    C = DummyExpression<dtype>(m, n);

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
    std::cout << "Temporary arena capacity: " << MemoryManager::temporary_arena().capacity()
              << " bytes, peak use: " << MemoryManager::temporary_arena().high_water_mark()
              << " bytes" << std::endl;
    
    

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