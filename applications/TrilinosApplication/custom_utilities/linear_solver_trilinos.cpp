// Project includes
#include "custom_utilities/linear_solver_trilinos.h"
#include "spaces/default_spaces.h"


namespace Kratos {


template class LinearSolver<
    TrilinosSpace<Epetra_FECrsMatrix,Epetra_FEVector>,
    TDefaultDenseSpace<double>
>;

template class TrilinosDofUpdater<
    TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>
>;


} // namespace Kratos
