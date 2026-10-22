#ifndef MATL_PRINT_SPECIALIZATIONS_HPP
#define MATL_PRINT_SPECIALIZATIONS_HPP

namespace matl 
{
   /***********************************************************************************************
    * void print(MAO const&) 
    * 
    * Prints a matrix, or the value of an expression, to the standard output stream, followed by a
    * newline (see operator<<).
    *
    * Parameters:
    *     mao         Matrix or expression to print
    ***********************************************************************************************/
    template<MatrixAlgebraObject MAO>
    void print(MAO const& mao)
    {std::cout << mao << std::endl;}
}
#endif
