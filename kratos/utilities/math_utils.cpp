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
    // The LU factorization works on local ublas dynamic copies so the same
    // code serves the (possibly Eigen-backed) Matrix/Vector interface types.
    const SizeType size1 = A.size1();
    DenseMatrix<TDataType> a_ublas(size1, A.size2());
    for (SizeType i = 0; i < size1; ++i)
        for (SizeType j = 0; j < A.size2(); ++j)
            a_ublas(i,j) = A(i,j);
    DenseVector<TDataType> x(rB.size());
    for (SizeType i = 0; i < static_cast<SizeType>(rB.size()); ++i)
        x[i] = rB[i];
    typedef permutation_matrix<SizeType> pmatrix;
    pmatrix pm(size1);
    int singular = lu_factorize(a_ublas,pm);
    KRATOS_DEBUG_ERROR_IF(singular == 1) << "Matrix is singular: " << A << std::endl;
    lu_substitute(a_ublas, pm, x);
    if (static_cast<SizeType>(rX.size()) != size1) rX.resize(size1, false);
    for (SizeType i = 0; i < size1; ++i)
        rX[i] = x[i];
}

/// Explicit instantation
template class MathUtils<double>;

} /// namespace Kratos