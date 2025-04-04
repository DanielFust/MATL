#include <iostream>
#include "../include/matl.hpp"






int main()
{
    using matl::RowMatrix;
    using matl::index_t;
    using matl::TypeTraits;
    using matl::_GENERIC_MAT_TRANSPOSE_;
    using matl::promote;

    RowMatrix<double> mat1(3,3);
    RowMatrix<float>  mat2(3,3);
    RowMatrix<int>    mat3(3,3);
    RowMatrix<float>  mat4(3,3);

    mat1 = 1;
    mat2 = 2.0;
    mat3 = 3;
    
   /* 
    * Scalar multiplication and matrix addition/subtraction do not require
    * dynamic allocation.
    * 
    * Different data types (e.g. float, double precision, integer) are
    * interoperable. Each consituent matrix operation uses native type 
    * promotion.
    */
    //mat4 = 1-mat1-mat2*transpose(mat1)-1; // No allocations
    //mat1 = mat1 + mat2*mat3;
    mat1 = mat2 + mat2*mat3;

   /*
    * Matrix multiplications detect if the LHS matrix is part of the RHS
    * expression and allocate temporary accordingly. 
    */ 
    //mat1 = mat1*mat2; // One allocation
    //mat4 = mat2*mat3; // No allocations

    //mat1.println();
    mat1.println();

    for (index_t i=1; i<=mat1.rows(); ++i)
        for (index_t j=1; j<=mat1.cols(); ++j)
            mat1(i,j) = i;
    mat1.println();
    mat1 = transpose(mat1);
    mat1.println();     
    //mat4.println();

    

    //TypeTraits< _GENERIC_MAT_TRANSPOSE_< RowMatrix<double> > >::dtype val;
    //TypeTraits< matl::_GENERIC_MAT_PROD_< RowMatrix<double>, RowMatrix<float> > >::dtype val;
    //_GENERIC_MAT_TRANSPOSE_< RowMatrix<double> >::dtype val;
    //promote<TypeTraits<RowMatrix<float> >::dtype,  TypeTraits<_GENERIC_MAT_TRANSPOSE_< RowMatrix<double> > >::dtype>::type val;

    float scalar = dot(mat1,transpose(mat1));   
    std::cout << "dot product value: " << scalar << std::endl;

    return 0;
}