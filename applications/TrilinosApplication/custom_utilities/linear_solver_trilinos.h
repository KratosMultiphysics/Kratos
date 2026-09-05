// Project includes
#include "linear_solvers/linear_solver.h"
#include "trilinos_space.h"
#include "spaces/ublas_space.h"
#include "spaces/default_spaces.h"

namespace Kratos {


KRATOS_API_EXTERN template class KRATOS_API(TRILINOS_APPLICATION) LinearSolver<
    TrilinosSpace<Epetra_FECrsMatrix,Epetra_FEVector>,
    TDefaultDenseSpace<double>
>;


} // namespace Kratos
