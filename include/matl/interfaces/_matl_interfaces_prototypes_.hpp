/**************************************************************************************************
 * File: _matl_interfaces_prototypes_.hpp
 *
 * Description: Base classes of the MATL matrix and expression types.
 *
 * The base classes use the curiously recurring template pattern (CRTP): a matrix type MAT derives
 * from MatrixBase<MAT>, and an expression type EXPR from ExpressionBase<EXPR>. They mark a type as
 * a matrix or an expression (see the MatrixType and ExpressionType concepts) and give generic code
 * access to the derived object. Generic functions are written as constrained templates, e.g.
 *
 *     template<MatrixType MAT> void f(MAT const& mat);
 *
 * so they work with the derived type directly, without virtual functions or forwarding.
 **************************************************************************************************/
#ifndef MATL_INTERFACES_PROTOTYPES_HPP
#define MATL_INTERFACES_PROTOTYPES_HPP

namespace matl 
{
   /***********************************************************************************************
    * class MatrixBase<MAT>
    * 
    * Base class of the MATL matrix types (e.g. RowMajorMatrix, FixedSizeRowMajorMatrix).
    * 
    * Template Parameters:
    *     MAT       Matrix type
    **********************************************************************************************/ 
    template<typename MAT>
    class MatrixBase
    {
      public:
        // Derived matrix
        MAT&       derived()       {return static_cast<MAT&>(*this);}
        MAT const& derived() const {return static_cast<MAT const&>(*this);}
    };


   /***********************************************************************************************
    * class ExpressionBase<EXPR>
    * 
    * Base class of the MATL expression types (see expression-types) and of views (MatrixView).
    * Every expression type provides:
    * 
    *     rows(), cols(), size()                   Dimensions of the result
    *     unsafe_evaluation_target(target)         'true' if evaluating directly into the target
    *                                              (a matrix or view) would overwrite elements
    *                                              still to be read
    *     reads(footprint)                         'true' if the expression may read elements in
    *                                              the given memory (see detail::Footprint)
    *     eval(target)                             Evaluate the expression into the target
    * 
    * Template Parameters:
    *     EXPR      Expression type
    **********************************************************************************************/ 
    template<typename EXPR>
    class ExpressionBase
    {
      public:
        // Derived expression
        EXPR const& derived() const {return static_cast<EXPR const&>(*this);}
    };
}
#endif
