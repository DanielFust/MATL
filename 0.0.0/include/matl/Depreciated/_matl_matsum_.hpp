#ifndef _MATL_MATSUM__HPP
#define _MATL_MATSUM__HPP
/**************************************************************************************************
 * Written by Daniel Fust
 *
 * This file contains template specializations of the MATSUM() function used to build matrix 
 * summation expression templates.
 **************************************************************************************************/
#include "_matl_forward_declr_.hpp"
#include "_matl_type_traits_.hpp"
#include "_matl_generic_mat_sum_.hpp"
#include "_matl_rowmat_sum_.hpp"

namespace matl
{


   /***********************************************************************************************
    *                            RowMatrix Summation Mappings
    ***********************************************************************************************/   

   /***********************************************************************************************
    * _ROWMAT_SUM_ MATSUM(RowMatrix const&, RowMatrix const&)
    * 
    * MATSUM template specialization for RowMatrix inputs
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
        MATSUM(RowMatrix<T1> const& mat1, RowMatrix<T2> const& mat2)
    {return _ROWMAT_SUM_<RowMatrix<T1>, RowMatrix<T2> > (mat1, mat2);}





   /**********************************************************************************************
    * _ROWMAT_SUM_ MATSUM(_ROWMAT_ATOMIC_EXPR_ const&, RowMatrix const&)
    * 
    * MATSUM template specialization for RowMatrix inputs
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
        MATSUM(_ROWMAT_ATOMIC_EXPR_<RMAE>  const& expr, RowMatrix<T> const& mat)
    {return _ROWMAT_SUM_< _ROWMAT_ATOMIC_EXPR_<RMAE>, RowMatrix<T> > (expr, mat);} 




   /**********************************************************************************************
    * _ROWMAT_SUM_ MATSUM(RowMatrix const&, _ROWMAT_ATOMIC_EXPR_ const&)
    * 
    * MATSUM template specialization for RowMatrix inputs
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
        MATSUM(RowMatrix<T> const& mat, _ROWMAT_ATOMIC_EXPR_<RMAE>  const& expr)
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
        MATSUM(_ROWMAT_ATOMIC_EXPR_<RMAE1>  const& expr1, _ROWMAT_ATOMIC_EXPR_<RMAE2>const& expr2)
    {
        return _ROWMAT_SUM_< _ROWMAT_ATOMIC_EXPR_<RMAE1>, 
                             _ROWMAT_ATOMIC_EXPR_<RMAE2> > (expr1, expr2);
    }   
    

    











   /***********************************************************************************************
    * _ROWMAT_SCALAR_SUM_ MATSUM(_ROWMAT_ATOMIC_EXPR_ const&, T);
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
    _ROWMAT_SCALAR_SUM_<RMAE,T> MATSUM(_ROWMAT_ATOMIC_EXPR_<RMAE> const& expr, T val)  
    {return _ROWMAT_SCALAR_SUM_<RMAE,T>(downcast(expr),val);} 




   /***********************************************************************************************
    * _ROWMAT_SCALAR_SUM_ MATSUM(T, _ROWMAT_ATOMIC_EXPR_ const&);
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
    _ROWMAT_SCALAR_SUM_<RMAE,T> MATSUM(T val,_ROWMAT_ATOMIC_EXPR_<RMAE> const& expr)  
    {return _ROWMAT_SCALAR_SUM_<RMAE,T>(downcast(expr),val);} 


  
   /***********************************************************************************************
    * _ROWMAT_SCALAR_SUM_ MATSUM(RowMatrix<T1> const&, T2);
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
    _ROWMAT_SCALAR_SUM_<RowMatrix<T1>,T2> MATSUM(RowMatrix<T1> const& mat, T2 val)  
    {return _ROWMAT_SCALAR_SUM_<RowMatrix<T1>,T2> (mat,val);}




   /***********************************************************************************************
    * _ROWMAT_SCALAR_SUM_ MATSUM(T2, RowMatrix<T1> const&);
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
    _ROWMAT_SCALAR_SUM_<RowMatrix<T1>,T2> MATSUM( T2 val, RowMatrix<T1> const& mat)  
    {return _ROWMAT_SCALAR_SUM_<RowMatrix<T1>,T2>(mat,val);} 



    

}
#endif