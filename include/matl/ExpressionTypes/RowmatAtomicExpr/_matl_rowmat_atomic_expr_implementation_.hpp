//#ifndef _MATL_ROWMAT_ATOMIC_EXPR_IMPLEMENTATION__HPP
//#define _MATL_ROWMAT_ATOMIC_EXPR_IMPLEMENTATION__HPP


//#include "_matl_forward_declarations_.hpp"
//#include "_matl_type_traits_.hpp"

namespace matl
{ 

   /***********************************************************************************************
    * void eval(Matrix& target) const
    * 
    * Evaluate expression to target matrix. Downcasts to derived type.
    * 
    *  Template Parameters
    *      EXPR_T           Atomic expression type
    *      MAT_T            Matrix type
    * 
    * Parameters
    *     Matrix target     Target matrix  
    ***********************************************************************************************/
    template<typename EXPR_T>
    template<typename MAT_T>
    void _ROWMAT_ATOMIC_EXPR_<EXPR_T>::eval(Matrix<MAT_T>& target) const
    {
        static_cast<EXPR_T const&>(*this).eval(static_cast<MAT_T&>(target));
        return;
    }


   /***********************************************************************************************
    * index_t rows() const
    * 
    * Number of rows in expression. Downcasts to derived type.
    * 
    *  Template Parameters:
    *      EXPR_T           Atomic expression type
    * 
    * Parameters:
    ***********************************************************************************************/
    template<typename EXPR_T>
    index_t _ROWMAT_ATOMIC_EXPR_<EXPR_T>::rows() const
    {return static_cast<EXPR_T const&>(*this).rows();}
    

   /***********************************************************************************************
    * index_t cols() const
    * 
    * Number of columns in expression. Downcasts to derived type.
    * 
    *  Template Parameters:
    *      EXPR_T           Atomic expression type
    * 
    * Parameters:
    ***********************************************************************************************/
    template<typename EXPR_T>
    index_t _ROWMAT_ATOMIC_EXPR_<EXPR_T>::cols() const
    {return static_cast<EXPR_T const&>(*this).cols();}
    


   /***********************************************************************************************
    * index_t size() const
    * 
    * Number of elements in expression. Downcasts to derived type.
    * 
    *  Template Parameters:
    *      EXPR_T           Atomic expression type
    * 
    * Parameters:
    ***********************************************************************************************/ 
    template<typename EXPR_T>
    index_t _ROWMAT_ATOMIC_EXPR_<EXPR_T>::size() const
    {return static_cast<EXPR_T const&>(*this).size();}

 
}


//#endif