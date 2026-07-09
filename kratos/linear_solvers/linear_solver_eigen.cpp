//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// Project includes
#include "linear_solvers/linear_solver_eigen.h"

// The Eigen-real LinearSolver instantiations exist only under the eigen
// backend: the two linear-algebra backends are mutually exclusive.
#ifdef KRATOS_USE_EIGEN_BACKEND

namespace Kratos {


template class LinearSolver<
    TEigenSparseSpace<double>,
    TUblasDenseSpace<double>
>;


template class LinearSolver<
    TEigenSparseSpace<float>,
    TUblasDenseSpace<double>
>;


} // namespace Kratos

#endif // KRATOS_USE_EIGEN_BACKEND
