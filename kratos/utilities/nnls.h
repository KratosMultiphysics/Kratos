//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Manuel Meßmer
//                   Ricky Aristio
//

#ifndef KRATOS_NNLS_INCLUDE_H
#define KRATOS_NNLS_INCLUDE_H

//// External includes
#include <boost/numeric/ublas/matrix.hpp>

namespace Kratos {
namespace nnls {

typedef boost::numeric::ublas::matrix<double> MatrixType;
typedef boost::numeric::ublas::vector<double> VectorType;

// Wrapper for nnls solver
double nnls(MatrixType& A, const VectorType& b, VectorType& x);

} // End Namespace nnls
} // End namespace Kratos

#endif // KRATOS_NNLS_INCLUDE_H