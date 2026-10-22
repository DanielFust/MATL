/**************************************************************************************************
 * File: _matl_elementwise_implementation_.hpp
 *
 * Description: Implementation of the elementwise expression template.
 **************************************************************************************************/
#ifndef MATL_ELEMENTWISE_IMPLEMENTATION_HPP
#define MATL_ELEMENTWISE_IMPLEMENTATION_HPP
namespace matl
{
    namespace expr
    {
       /*******************************************************************************************
        * Constructor
        *
        * Stores the operands: expressions and scalars are copied, matrices are referenced.
        *
        * Parameters:
        *     lhs          Left operand
        *     rhs          Right operand
        *******************************************************************************************/
        template<typename OP, typename LHS, typename RHS>
        Elementwise<OP,LHS,RHS>::Elementwise(LHS const& lhs, RHS const& rhs) : lhs_(lhs), rhs_(rhs)
        {}

       /*******************************************************************************************
        * LHS const& lhs() const
        *
        * Returns:
        *     Left operand
        *******************************************************************************************/
        template<typename OP, typename LHS, typename RHS>
        LHS const& Elementwise<OP,LHS,RHS>::lhs() const
        {return this->lhs_;}

       /*******************************************************************************************
        * RHS const& rhs() const
        *
        * Returns:
        *     Right operand
        *******************************************************************************************/
        template<typename OP, typename LHS, typename RHS>
        RHS const& Elementwise<OP,LHS,RHS>::rhs() const
        {return this->rhs_;}

       /*******************************************************************************************
        * index_t rows() const
        *
        * Returns:
        *     index_t      Number of rows in resultant matrix
        *******************************************************************************************/
        template<typename OP, typename LHS, typename RHS>
        index_t Elementwise<OP,LHS,RHS>::rows() const
        {
            if constexpr (Scalar<LHS>) return this->rhs_.rows();
            else                       return this->lhs_.rows();
        }

       /*******************************************************************************************
        * index_t cols() const
        *
        * Returns:
        *     index_t      Number of columns in resultant matrix
        *******************************************************************************************/
        template<typename OP, typename LHS, typename RHS>
        index_t Elementwise<OP,LHS,RHS>::cols() const
        {
            if constexpr (Scalar<LHS>) return this->rhs_.cols();
            else                       return this->lhs_.cols();
        }

       /*******************************************************************************************
        * index_t size() const
        *
        * Returns:
        *     index_t      Number of elements in resultant matrix
        *******************************************************************************************/
        template<typename OP, typename LHS, typename RHS>
        index_t Elementwise<OP,LHS,RHS>::size() const
        {
            if constexpr (Scalar<LHS>) return this->rhs_.size();
            else                       return this->lhs_.size();
        }

       /*******************************************************************************************
        * void check_dimensions() const
        *
        * Checks the operands (recursively), then that matrix operands have the same dimensions.
        * Reports a matl::dimension_error otherwise.
        *******************************************************************************************/
        template<typename OP, typename LHS, typename RHS>
        void Elementwise<OP,LHS,RHS>::check_dimensions() const
        {
            detail::check_dimensions(this->lhs_);
            detail::check_dimensions(this->rhs_);

            if constexpr (!Scalar<LHS> && !Scalar<RHS>)
                if ((this->lhs_.rows() != this->rhs_.rows()) || (this->lhs_.cols() != this->rhs_.cols()))
                    detail::incompatible_dimensions("elementwise operation",
                                            this->lhs_.rows(), this->lhs_.cols(),
                                            this->rhs_.rows(), this->rhs_.cols());
        }

       /*******************************************************************************************
        * bool unsafe_evaluation_target(detail::Footprint const&)
        *
        * Check if evaluation of the expression directly to the target matrix would corrupt
        * memory. Elementwise evaluation is unsafe if the target is read through a transpose
        * (element (i,j) of the result would read an element of the target that has already been
        * overwritten), or if an operand is unsafe on its own.
        *
        * Parameters:
        *     target       Footprint of the target matrix or view
        * Returns:
        *     'true' if evaluation directly to the target matrix would corrupt memory
        *     'false' otherwise
        *******************************************************************************************/
        template<typename OP, typename LHS, typename RHS>
        bool Elementwise<OP,LHS,RHS>::unsafe_evaluation_target(detail::Footprint const& target) const
        {
            return detail::unsafe_elementwise_operand(target, this->lhs_) ||
                   detail::unsafe_elementwise_operand(target, this->rhs_);
        }

       /*******************************************************************************************
        * bool reads(detail::Footprint const&)
        *
        * Check if the expression may read elements of the target (see detail::Footprint).
        *
        * Parameters:
        *     target       Footprint of the target (where its elements are in memory)
        * Returns:
        *     'true' if the expression references the target matrix, 'false' otherwise
        *******************************************************************************************/
        template<typename OP, typename LHS, typename RHS>
        bool Elementwise<OP,LHS,RHS>::reads(detail::Footprint const& target) const
        {return detail::references(target, this->lhs_) || detail::references(target, this->rhs_);}

       /*******************************************************************************************
        * void eval(MAT& target) const
        *
        * Evaluate the expression to the target matrix.
        *
        * Parameters:
        *     target       Target matrix (must be a safe evaluation target with matching dimensions)
        *******************************************************************************************/
        template<typename OP, typename LHS, typename RHS>
        template<WritableMatrix MAT>
        void Elementwise<OP,LHS,RHS>::eval(MAT& target) const
        {detail::assign_elementwise(target, *this);}
    }

    namespace detail
    {
       /*******************************************************************************************
        * struct evaluator<expr::Elementwise>
        *
        * Evaluates elements of the expression from the evaluators of its operands.
        *******************************************************************************************/
        template<typename OP, typename LHS, typename RHS>
        struct evaluator< expr::Elementwise<OP,LHS,RHS> >
        {
            typedef expr::Elementwise<OP,LHS,RHS>  expression_type;
            typedef typename expression_type::dtype dtype;

            evaluator<LHS> lhs_;
            evaluator<RHS> rhs_;

            explicit evaluator(expression_type const& expr) : lhs_(expr.lhs()), rhs_(expr.rhs()) {}

            dtype operator() (index_t i, index_t j) const {return OP::apply(lhs_(i,j), rhs_(i,j));}
            dtype operator() (index_t i)            const {return OP::apply(lhs_(i), rhs_(i));}
        };
    }
}
#endif
