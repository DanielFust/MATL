#ifndef _MATL_OPERATOR_PLUS__HPP
#define _MATL_OPERATOR_PLUS__HPP
/**************************************************************************************************
 * Written by Daniel Fust
 *
 * This file contains implementation of the operator(+) for the BLAS library
 **************************************************************************************************/
#include "_matl_forward_declr_.hpp"
#include "_matl_type_traits_.hpp"
#include "_matl_matsum_.hpp"
#include "_matl_rowmat_scalar_sum_.hpp"

namespace matl
{

   //#if 0

   /***********************************************************************************************
    * _ATOMIC_EXPR_ operator+ (Matrix const&, _MATEXPRESSION_ const&);
    * 
    * General summation of matrix and matrix expression
    * 
    * Template Parameters:
    *     EXPR_T        (implicit) Expression type
    *     MAT_T         (implicit) matrix type
    * 
    * Parameters:
    *     expr           Matrix expression
    *     mat            Matrix
    * 
    * Returns:
    *     _ATOMIC_EXPR_ derived type atomic matrix expression superclass
    ***********************************************************************************************/ 
    template<typename EXPR_T, typename MAT_T>
    auto operator+ (_MATEXPRESSION_<EXPR_T> const& expr,
                    Matrix<MAT_T>           const& mat)  
    {return MATSUM(downcast(expr),downcast(mat));}






   /***********************************************************************************************
    * _ATOMIC_EXPR_ operator+ (_MATEXPRESSION_ const&, Matrix const&,);
    * 
    * General summation of matrix expression and matrix
    * 
    * Template Parameters:
    *     MAT_T         (implicit) matrix type
    *     EXPR_T        (implicit) Expression type
    * 
    * Parameters:
    *     mat            Matrix
    *     expr           Matrix expression
    * 
    * Returns:
    *     _ATOMIC_EXPR_ derived type atomic matrix expression superclass
    ***********************************************************************************************/ 
    template<typename MAT_T, typename EXPR_T>
    auto operator+ (Matrix<MAT_T>           const& mat,
                    _MATEXPRESSION_<EXPR_T> const& expr)
    {return MATSUM(downcast(mat),downcast(expr));}







   /***********************************************************************************************
    * _ATOMIC_EXPR_ operator+ (Matrix const&, Matrix const&,);
    * 
    * General summation of matrices
    * 
    * Template Parameters:
    *     MAT_T1         (implicit) first matrix type
    *     MAT_T2         (implicit) second matrix type
    * 
    * Parameters:
    *     mat1            First matrix
    *     mat2            Second matrix
    * 
    * Returns:
    *     _ATOMIC_EXPR_ derived type atomic matrix expression superclass
    ***********************************************************************************************/ 
    template<typename MAT_T1, typename MAT_T2>
    auto operator+ (Matrix<MAT_T1>           const& mat1,
                    Matrix<MAT_T2>           const& mat2)
    {return MATSUM(downcast(mat1),downcast(mat2));}  

 





   /***********************************************************************************************
    * _ATOMIC_EXPR_ operator+ (_MATEXPRESSION_ const&, const&, _MATEXPRESSION_ const&);
    * 
    * General summation of expressions
    * 
    * Template Parameters:
    *     EXPR_T1       (implicit) First expression type
    *     EXPR_T2       (implicit) Second expression type
    * 
    * Parameters:
    *     expr1          First matrix expression
    *     expr2          Second matrix expression
    * 
    * Returns:
    *     _ATOMIC_EXPR_ derived type atomic matrix expression superclass
    ***********************************************************************************************/ 
    template<typename EXPR_T1, typename EXPR_T2>
    auto operator+ (_MATEXPRESSION_<EXPR_T1> const& expr1,
                    _MATEXPRESSION_<EXPR_T2> const& expr2) 
    {return MATSUM(downcast(expr1), downcast(expr2));} 










   /***********************************************************************************************
    * _ATOMIC_EXPR_ operator+ (_MATEXPRESSION_ const&, T);
    * 
    * General summation of expression and scalar
    * 
    * Template Parameters:
    *     EXPR_T       (implicit) Expression type
    *     T            (implicit) Scalar data type
    * 
    * Parameters:
    *     expr           Matrix expression
    *     val            Scalar
    * 
    * Returns:
    *     _ATOMIC_EXPR_ derived type atomic matrix expression superclass
    ***********************************************************************************************/ 
    template<typename EXPR_T, class T, enable_if_scalar<T> = true>
    auto operator+ (_MATEXPRESSION_<EXPR_T> const& expr, T val) 
    {return MATSUM(downcast(expr), val);}



   /***********************************************************************************************
    * _ATOMIC_EXPR_ operator+ (T, _MATEXPRESSION_ const&);
    * 
    * General summation of expression and scalar
    * 
    * Template Parameters:
    *     T            (implicit) Scalar data type
    *     EXPR_T       (implicit) Expression type
    * 
    * Parameters:
    *     val            Scalar
    *     expr           Matrix expression
    * 
    * Returns:
    *     _ATOMIC_EXPR_ derived type atomic matrix expression superclass
    ***********************************************************************************************/ 
    template<typename EXPR_T, class T, enable_if_scalar<T> = true>
    auto operator+ (T val, _MATEXPRESSION_<EXPR_T> const& expr) 
    {return MATSUM(val, downcast(expr));}   


   /***********************************************************************************************
    * _ATOMIC_EXPR_ operator+ (Matrix const&, T);
    * 
    * General summation of Matrix and scalar
    * 
    * Template Parameters:
    *     MAT_T        (implicit) Matrix type
    *     T            (implicit) Scalar data type
    * 
    * Parameters:
    *     mat           Matrix expression
    *     val           Scalar
    * 
    * Returns:
    *     _ATOMIC_EXPR_ derived type atomic matrix expression superclass
    ***********************************************************************************************/ 
    template<typename MAT_T, class T, enable_if_scalar<T> = true>
    auto operator+ (Matrix<MAT_T> const& mat, T val) 
    {return MATSUM(downcast(mat), val);} 





   /***********************************************************************************************
    * _ATOMIC_EXPR_ operator+ (T, Matrix const&);
    * 
    * General summation of Matrix and scalar
    * 
    * Template Parameters:
    *     MAT_T        (implicit) Matrix type
    *     T            (implicit) Scalar data type
    * 
    * Parameters:
    *     val           Scalar
    *     mat           Matrix expression
    * 
    * Returns:
    *     _ATOMIC_EXPR_ derived type atomic matrix expression superclass
    ***********************************************************************************************/ 
    template<typename MAT_T, class T, enable_if_scalar<T> = true>
    auto operator+ (T val, Matrix<MAT_T> const& mat) 
    {return MATSUM(val, downcast(mat));} 







    //#endif


    #if 0







   /***********************************************************************************************
    *                            RowMatrix Summation Mappings
    ***********************************************************************************************/   



   /***********************************************************************************************
    * _ROWMAT_SCALAR_SUM_ operator+ (RowMatrix const&, T);
    * 
    * General summation of Matrix and scalar
    * 
    * Template Parameters:
    *     T1           (implicit) RowMatrix data type
    *     T2           (implicit) Scalar data type
    * 
    * Parameters:
    *     mat           RowMatrix 
    *     val           Scalar
    * 
    * Returns:
    *     _ROWMAT_SCALAR_EXPR_ derived type atomic matrix expression superclass
    ***********************************************************************************************/ 
    template<class T1, class T2>
    auto operator+ (RowMatrix<T1> const& mat, T2 val) 
    {return MATSUM(downcast(mat), val);} 





   /***********************************************************************************************
    * _ROWMAT_SCALAR_SUM_ operator+ (T, RowMatrix const&);
    * 
    * General summation of RowMatrix and scalar
    * 
    * Template Parameters:
    *     T1           (implicit) RowMatrix data type
    *     T2           (implicit) Scalar data type
    * 
    * Parameters:
    *     val           Scalar
    *     mat           RowMatrix
    * 
    * Returns:
    *     _ATOMIC_EXPR_ derived type atomic matrix expression superclass
    ***********************************************************************************************/ 
    template<class T1, class T2>
    auto operator+ (T2 val, RowMatrix<T1> const& mat) 
    {return MATSUM(val, downcast(mat));} 



    #endif


















    #if 0




   /***********************************************************************************************
    *                            RowMatrix Summation Mappings
    ***********************************************************************************************/   

   /***********************************************************************************************
    * _ROWMAT_SUM_ operator+(RowMatrix const&, RowMatrix const&)
    * 
    * Addition operator for RowMatrix inputs
    * 
    * Template Parameters:
    *     T1 (implicit)      First RowMatrix data type
    *     T2 (implicit)      Second RowMatrix data type
    * 
    * Parameters:
    *     mat1               First RowMatrix
    *     mat2               Second RowMatrix
    * 
    * Returns:
    *     _ROWMAT_SUM_       Row matrix summation expression template
    **********************************************************************************************/ 
    template<class T1, class T2>
    _ROWMAT_SUM_<RowMatrix<T1>, RowMatrix<T2> > 
         operator+(RowMatrix<T1> const& mat1, RowMatrix<T2> const& mat2)
    {return _ROWMAT_SUM_<RowMatrix<T1>, RowMatrix<T2> > (mat1, mat2);}





   /**********************************************************************************************
    * _ROWMAT_SUM_  operator+(_ROWMAT_ATOMIC_EXPR_ const&, RowMatrix const&)
    * 
    * Addition operator for for Row Matrix Expression - RowMatrix inputs
    * 
    * Template Parameters:
    *     RMAE (implicit)     Row matrix atomic expression type data type
    *     T    (implicit)     RowMatrix data type
    * 
    * Parameters:
    *     expr               Row matrix atomic expression
    *     mat                RowMatrix
    * 
    * Returns:
    *     _ROWMAT_SUM_       Row matrix summation expression template
    **********************************************************************************************/ 
    template<class RMAE, class T>
    _ROWMAT_SUM_< _ROWMAT_ATOMIC_EXPR_<RMAE>, RowMatrix<T> > 
         operator+(_ROWMAT_ATOMIC_EXPR_<RMAE>  const& expr, RowMatrix<T> const& mat)
    {return _ROWMAT_SUM_< RMAE, RowMatrix<T> > (downcast(expr), mat);} 




   /**********************************************************************************************
    * _ROWMAT_SUM_  operator+(RowMatrix const&, _ROWMAT_ATOMIC_EXPR_ const&)
    * 
    * Addition operator for RowMatrix Row Matrix expression inputs
    * 
    * Template Parameters:
    *     T    (implicit)     RowMatrix data type
    *     RMAE (implicit)     Row matrix atomic expression type data type
    * 
    * Parameters:
    *     mat                RowMatrix
    *     expr               Row matrix atomic expression
    * 
    * Returns:
    *     _ROWMAT_SUM_       Row matrix summation expression template
    **********************************************************************************************/ 
    template<class T, class RMAE>
    _ROWMAT_SUM_< RowMatrix<T>, _ROWMAT_ATOMIC_EXPR_<RMAE> > 
         operator+(RowMatrix<T> const& mat, _ROWMAT_ATOMIC_EXPR_<RMAE>  const& expr)
    {return _ROWMAT_SUM_< RowMatrix<T>, _ROWMAT_ATOMIC_EXPR_<RMAE> > (mat, expr);}  




   /**********************************************************************************************
    * _ROWMAT_SUM_ MATSUM(_ROWMAT_ATOMIC_EXPR_ const&, _ROWMAT_ATOMIC_EXPR_ const&)
    * 
    * MATSUM template specialization for RowMatrix inputs
    * 
    * Template Parameters:
    *     RMAE1 (implicit)     First row matrix atomic expression type data type
    *     RMAE2 (implicit)     Second row matrix atomic expression type data type
    * 
    * Parameters:
    *     expr1                Row matrix atomic expression
    *     expr2                Row matrix atomic expression
    * 
    * Returns:
    *     _ROWMAT_SUM_         Row matrix summation expression template
    **********************************************************************************************/ 
    template<class RMAE1, class RMAE2>
    _ROWMAT_SUM_< _ROWMAT_ATOMIC_EXPR_<RMAE1>, _ROWMAT_ATOMIC_EXPR_<RMAE2> > 
         operator+(_ROWMAT_ATOMIC_EXPR_<RMAE1>  const& expr1, _ROWMAT_ATOMIC_EXPR_<RMAE2>const& expr2)
    {
        return _ROWMAT_SUM_< _ROWMAT_ATOMIC_EXPR_<RMAE1>, 
                             _ROWMAT_ATOMIC_EXPR_<RMAE2> > (expr1, expr2);
    }   
    

    






   /***********************************************************************************************
    * _ROWMAT_SCALAR_SUM_  operator+(_ROWMAT_ATOMIC_EXPR_ const&, T);
    * 
    * Addition of Row Matrix Atomic Expression and scalar
    * 
    * Template Parameters:
    *     RMAE        (implicit) Row Matrix Atomic Expression type
    *     T           (implicit) Scalar data type
    * 
    * Parameters:
    *     expr         Matrix expression
    *     val          Scalar
    * 
    * Returns:
    *     _ROWMAT_SCALAR_SUM_ Row matrix scalar product
    ***********************************************************************************************/ 
    template<typename RMAE, class T>
    _ROWMAT_SCALAR_SUM_<RMAE,T>  operator+(_ROWMAT_ATOMIC_EXPR_<RMAE> const& expr, T val)  
    {return _ROWMAT_SCALAR_SUM_<RMAE,T>(downcast(expr),val);} 




   /***********************************************************************************************
    * _ROWMAT_SCALAR_SUM_  operator+(T, _ROWMAT_ATOMIC_EXPR_ const&);
    * 
    * Addition of Row Matrix Atomic Expression and scalar
    * 
    * Template Parameters:
    *     RMAE        (implicit) Row Matrix Atomic Expression type
    *     T           (implicit) Scalar data type
    * 
    * Parameters:
    *     expr         Matrix expression
    *     val          Scalar
    * 
    * Returns:
    *     _ROWMAT_SCALAR_SUM_ Row matrix scalar product
    ***********************************************************************************************/ 
    template<typename RMAE, class T>
    _ROWMAT_SCALAR_SUM_<RMAE,T>  operator+(T val,_ROWMAT_ATOMIC_EXPR_<RMAE> const& expr)  
    {return _ROWMAT_SCALAR_SUM_<RMAE,T>(downcast(expr),val);} 


  
   /***********************************************************************************************
    * _ROWMAT_SCALAR_SUM_  operator+(RowMatrix<T1> const&, T2);
    * 
    * Addition of Row Matrix Atomic Expression and scalar
    * 
    * Template Parameters:
    *     T1          (implicit) Row Matrix data type
    *     T2          (implicit) Scalar data type
    * 
    * Parameters:
    *     mat          Row Matrix
    *     val          Scalar
    * 
    * Returns:
    *     _ROWMAT_SCALAR_SUM_ Row matrix scalar product
    ***********************************************************************************************/ 
    template<typename T1, class T2>
    _ROWMAT_SCALAR_SUM_<RowMatrix<T1>,T2>  operator+(RowMatrix<T1> const& mat, T2 val)  
    {return _ROWMAT_SCALAR_SUM_<RowMatrix<T1>,T2> (mat,val);}




   /***********************************************************************************************
    * _ROWMAT_SCALAR_SUM_  operator+(T2, RowMatrix<T1> const&);
    * 
    * Addition of Row Matrix Atomic Expression and scalar
    * 
    * Template Parameters:
    *     T1          (implicit) Row Matrix data type
    *     T2          (implicit) Scalar data type
    * 
    * Parameters:
    *     val          Scalar
    *     mat          Row Matrix
    * 
    * Returns:
    *     _ROWMAT_SCALAR_SUM_ Row matrix scalar product
    ***********************************************************************************************/ 
    template<typename T1, class T2>
    _ROWMAT_SCALAR_SUM_<RowMatrix<T1>,T2>  operator+( T2 val, RowMatrix<T1> const& mat)  
    {return _ROWMAT_SCALAR_SUM_<RowMatrix<T1>,T2>(mat,val);} 






    #endif














    




}
#endif