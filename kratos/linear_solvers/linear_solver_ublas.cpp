// Project includes
#include "linear_solvers/linear_solver_ublas.h"


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

template class Reorderer<
    TUblasSparseSpace<double>,
    TUblasDenseSpace<double>
>;

template class Reorderer<
    TUblasSparseSpace<float>,
    TUblasDenseSpace<double>
>;

template class DofUpdater<
    TUblasSparseSpace<double>
>;

template class DofUpdater<
    TUblasSparseSpace<float>
>;

} // namespace Kratos
