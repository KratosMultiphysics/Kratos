//    |  /           | 
//    ' /   __| _` | __|  _ \   __| 
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/ 
//                   Multi-Physics  
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//

#if !defined(KRATOS_AMATRIX_INTERFACE_H_INCLUDED )
#define  KRATOS_AMATRIX_INTERFACE_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes
#include "amatrix.h"

// Project includes
#include "includes/define.h"
#include "includes/dense_matrix.h"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

using Matrix = DenseMatrix<double,AMatrix::dynamic, AMatrix::dynamic>;

using Vector = DenseMatrix<double,AMatrix::dynamic, 1>;

template <typename T> T& noalias(T& TheMatrix){return TheMatrix.noalias();}

template <typename T> AMatrix::TransposeMatrix<T> trans(T& TheMatrix){return TheMatrix.transpose();}

using ZeroMatrix = AMatrix::ZeroMatrix<double>;

template <typename TExpression1Type, typename TExpression2Type,
    std::size_t TCategory1, std::size_t TCategory2>
AMatrix::MatrixProductExpression<TExpression1Type, TExpression2Type> prod(
    AMatrix::MatrixExpression<TExpression1Type, TCategory1> const& First,
    AMatrix::MatrixExpression<TExpression2Type, TCategory2> const& Second) {
    return AMatrix::MatrixProductExpression<TExpression1Type, TExpression2Type>(
        First.expression(), Second.expression());
}


///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{


///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

}  // namespace Kratos.

#endif // KRATOS_AMATRIX_INTERFACE_H_INCLUDED  defined 



