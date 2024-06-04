//   ____        ______  _____
//  / __ \      |  ____|/ ____|
// | |  | |_   _| |__  | (___   ___
// | |  | | | | |  __|  \___ \ / _ \'
// | |__| | |_| | |____ ____) | (_) |
//  \___\_\\__,_|______|_____/ \___/
//         Quadrature for Embedded Solids
//
//  License:    BSD 4-Clause License
//              See: https://github.com/manuelmessmer/QuESo/blob/main/LICENSE
//
//  Authors:    Manuel Messmer


#ifndef NNLS_INCLUDE_H
#define NNLS_INCLUDE_H

//// STL includes
#include <vector>

namespace queso {

///@name QuESo Classes
///@{


/**
 * @class  NNLS solver
 * @brief  Wrapper to Non-Negative-Least-Squares-Solver (see: external_libaries/nnls/nnls_impl.h).
 * @author Manuel Messmer
**/
class NNLS {
public:

    ///@}
    ///@name Type Definitions
    ///@{

    typedef std::vector<double> MatrixType;
    typedef std::vector<double> VectorType;

    ///@}
    ///@name Operations
    ///@{

    /// @brief  Solves ||Ax-b||_L2 as a non-negative least-squares problem. All entries in x are positive.
    ///         For better performance A and B are not copied. However, this means they are non-const!!
    /// @param A Matrix as vector form (Serialized matrix, column first)
    /// @param b RHS
    /// @param X Solution Vector
    /// @return Residuum (double)
    static double solve(MatrixType& A, VectorType& B, VectorType& X);

    ///@}

}; // End Class NNLS
///@}
} // End namespace queso

#endif // NNLS_INCLUDE_H