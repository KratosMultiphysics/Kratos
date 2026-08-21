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
