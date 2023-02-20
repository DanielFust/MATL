#include <iostream>
#include "../include/matl.hpp"




int main()
{
    using matl::RowMatrix;

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
    mat4 = 5*mat1*mat2*mat2; // No allocations

   /*
    * Matrix multiplications detect if the LHS matrix is part of the RHS
    * expression and allocate temporary accordingly. 
    */ 
    //mat1 = mat1*mat2; // One allocation
    //mat4 = mat2*mat3; // No allocations

    mat1.println();
    mat4.println();
    return 0;
}