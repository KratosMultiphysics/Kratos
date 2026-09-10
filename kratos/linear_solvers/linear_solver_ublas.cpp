// Project includes
#include "linear_solvers/linear_solver_ublas.h"

// The uBLAS-real LinearSolver instantiations exist only under the ublas
// backend: the two linear-algebra backends are mutually exclusive.
#ifndef KRATOS_USE_EIGEN_BACKEND

namespace Kratos {


template class LinearSolver<
    TUblasSparseSpace<double>,
    TUblasDenseSpace<double>
>;


template class LinearSolver<
    TUblasSparseSpace<float>,
    TUblasDenseSpace<double>
>;


} // namespace Kratos

#endif // !KRATOS_USE_EIGEN_BACKEND
