/**************************************************************************************************
 * File: _matl_unary_implementation_.hpp
 *
 * Description: Implementation of the unary elementwise expression template.
 **************************************************************************************************/
#ifndef MATL_UNARY_IMPLEMENTATION_HPP
#define MATL_UNARY_IMPLEMENTATION_HPP
namespace matl
{
    namespace expr
    {
       /*******************************************************************************************
        * Constructor
        *
        * Stores the operand: an expression is copied, a matrix is referenced.
        *
        * Parameters:
        *     mat          Operand (matrix or expression)
        *******************************************************************************************/
        template<typename OP, typename MAO>
        Unary<OP,MAO>::Unary(MAO const& mat) : mat_(mat)
        {}

       /*******************************************************************************************
        * MAO const& operand() const
        *
        * Returns:
        *     Operand
        *******************************************************************************************/
        template<typename OP, typename MAO>
        MAO const& Unary<OP,MAO>::operand() const
        {return this->mat_;}

       /*******************************************************************************************
        * index_t rows() const, cols() const, size() const
        *
        * Returns:
        *     Dimensions of the result (those of the operand)
        *******************************************************************************************/
        template<typename OP, typename MAO>
        index_t Unary<OP,MAO>::rows() const
        {return this->mat_.rows();}

        template<typename OP, typename MAO>
        index_t Unary<OP,MAO>::cols() const
        {return this->mat_.cols();}

        template<typename OP, typename MAO>
        index_t Unary<OP,MAO>::size() const
        {return this->mat_.size();}

       /*******************************************************************************************
        * void check_dimensions() const
        *
        * A unary operation fits any operand; checks the operand (recursively).
        *******************************************************************************************/
        template<typename OP, typename MAO>
        void Unary<OP,MAO>::check_dimensions() const
        {detail::check_dimensions(this->mat_);}

       /*******************************************************************************************
        * bool unsafe_evaluation_target(MAT const&)
        *
        * Element (i,j) of the result depends only on element (i,j) of the operand, so evaluation
        * into the target is unsafe only if the operand reads the target through a transpose or is
        * unsafe itself (see detail::unsafe_elementwise_operand).
        *
        * Parameters:
        *     target       Target matrix
        * Returns:
        *     'true' if evaluation directly to the target matrix would corrupt memory
        *******************************************************************************************/
        template<typename OP, typename MAO>
        template<WritableMatrix MAT>
        bool Unary<OP,MAO>::unsafe_evaluation_target(MAT const& target) const
        {return detail::unsafe_elementwise_operand(target, this->mat_);}

       /*******************************************************************************************
        * bool reads(detail::Footprint const&)
        *
        * Returns:
        *     'true' if the expression may read elements of the target (see detail::Footprint)
        *******************************************************************************************/
        template<typename OP, typename MAO>
        bool Unary<OP,MAO>::reads(detail::Footprint const& target) const
        {return detail::references(target, this->mat_);}

       /*******************************************************************************************
        * void eval(MAT& target) const
        *
        * Evaluate the expression to the target matrix.
        *
        * Parameters:
        *     target       Target matrix (must be a safe evaluation target with matching dimensions)
        *******************************************************************************************/
        template<typename OP, typename MAO>
        template<WritableMatrix MAT>
        void Unary<OP,MAO>::eval(MAT& target) const
        {detail::assign_elementwise(target, *this);}
    }

    namespace detail
    {
       /*******************************************************************************************
        * struct evaluator<expr::Unary>
        *
        * Evaluates elements of the expression from the evaluator of its operand.
        *******************************************************************************************/
        template<typename OP, typename MAO>
        struct evaluator< expr::Unary<OP,MAO> >
        {
            typedef expr::Unary<OP,MAO>             expression_type;
            typedef typename expression_type::dtype dtype;

            evaluator<MAO> operand_;

            explicit evaluator(expression_type const& expr) : operand_(expr.operand()) {}

            dtype operator() (index_t i, index_t j) const {return OP::apply(operand_(i,j));}
            dtype operator() (index_t i)            const {return OP::apply(operand_(i));}
        };
    }
}
#endif
