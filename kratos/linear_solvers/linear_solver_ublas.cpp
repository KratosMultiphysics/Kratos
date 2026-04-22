// Project includes
#include "linear_solvers/linear_solver_ublas.h"


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
