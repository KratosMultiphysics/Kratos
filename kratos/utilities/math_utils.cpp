//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//
//  Collaborators:   Vicente Mataix Ferrandiz
//                   Pablo Becker
//

// System includes

// External includes

// Project includes
#include "utilities/math_utils.h"

namespace Kratos
{

template<class TDataType>
void MathUtils<TDataType>::Solve(
    MatrixType A,
    VectorType& rX,
    const VectorType& rB
    )
{
    const SizeType size1 = A.size1();
    rX = rB;
    typedef permutation_matrix<SizeType> pmatrix;
    pmatrix pm(size1);
    int singular = lu_factorize(A,pm);
    KRATOS_DEBUG_ERROR_IF(singular == 1) << "Matrix is singular: " << A << std::endl;
    lu_substitute(A, pm, rX);
}

/// Explicit instantation
template class MathUtils<double>;

} /// namespace Kratos