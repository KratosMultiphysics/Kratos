// Project includes
#include "custom_utilities/linear_solver_trilinos.h"


namespace Kratos {


template class LinearSolver<
    TrilinosSpace<Epetra_FECrsMatrix,Epetra_FEVector>,
    UblasSpace<double, Matrix, Vector>
>;


} // namespace Kratos
