// The uBLAS-real LinearSolver instantiations exist only under the ublas
// backend: the two linear-algebra backends are mutually exclusive.
#ifndef KRATOS_USE_EIGEN_BACKEND

// Project includes
#include "linear_solvers/linear_solver.h"
#include "spaces/ublas_space.h"

namespace Kratos {


KRATOS_API_EXTERN template class KRATOS_API(KRATOS_CORE) LinearSolver<
    TUblasSparseSpace<double>,
    TUblasDenseSpace<double>
>;


KRATOS_API_EXTERN template class KRATOS_API(KRATOS_CORE) LinearSolver<
    TUblasSparseSpace<float>,
    TUblasDenseSpace<double>
>;


} // namespace Kratos

#endif // !KRATOS_USE_EIGEN_BACKEND
