#include <iostream>



class MyClass {};

template<typename MAT_T>
class Matrix {};

template<class T>
class RowMatrix : public Matrix<RowMatrix<T> > 
{
    T val_;

  public:  
    
    RowMatrix() {this->val_=0;}

    template<class T2>
    RowMatrix(T2 val);

    //RowMatrix(T val);
}; 

//template<class T>
//RowMatrix<T>::RowMatrix(T val){};

template<class T1>
template<class T2>
RowMatrix<T1>::RowMatrix(T2 val)
{
    std::cout << "Generic constructor called" << std::endl;
};



template<typename ...T> struct TypeTraits;

template<class T>
struct TypeTraits<RowMatrix<T> >
{
    typedef RowMatrix<T> downcast_type;
}; 

template<class T>
struct TypeTraits<T>
{
    typedef T downcast_type;
}; 




// Enable if MyType
template<typename T>
//using enable_if_my_type = std::enable_if_t<std::is_same<T, MyClass>::value >;
using enable_if_my_type = std::enable_if_t<std::is_same<T, MyClass>::value, bool>;


// Enable if Matrix (Doesn't work)
template<typename MAT_T>
using enable_if_matrix = std::enable_if_t<std::is_same<MAT_T, Matrix<typename TypeTraits<MAT_T>::downcast_type> >::value, bool>;

template<typename MAT_T>
using is_matrix = std::is_same< MAT_T, Matrix<typename TypeTraits<MAT_T>::downcast_type > >;











// Enable if MyType
template<typename T, enable_if_my_type<T> = true>
void my_func(T& obj)
{
    std::cout << "Yeah, that's my type" << std::endl;
}


template<typename MAT_T>
Matrix<MAT_T> my_func2(Matrix<MAT_T> const& mat)
{
    MAT_T mat2;
    return mat2;
}



int main()
{
    MyClass obj;
    int     a = 1;

    //my_func(obj);
    //my_func(a);

    std::cout << "is_matrix<RowMatrix<float> >::value = " << is_matrix<RowMatrix<float> >::value << std::endl;

    std::cout << "std::is_same<MyClass, MyClass>::value = " << std::is_same<MyClass, MyClass>::value << std::endl;


    RowMatrix<double> mat;
    auto mat2 = my_func2(mat);

    return 0;
}



