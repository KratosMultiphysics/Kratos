// Project includes
#include "linear_solvers/linear_solver.h"
#include "spaces/ublas_space.h"


namespace Kratos {

template class LinearSolver<
    TUblasSparseSpace<double>,
    TUblasDenseSpace<double>,
    Reorderer<TUblasSparseSpace<double>, TUblasDenseSpace<double>>
>;

template class LinearSolver<
    TUblasSparseSpace<float>,
    TUblasDenseSpace<double>,
    Reorderer<TUblasSparseSpace<float>, TUblasDenseSpace<double>>
>;

} // namespace Kratos
