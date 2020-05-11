// System includes

#ifdef KRATOS_USE_AMATRIX
#include "boost/numeric/ublas/matrix.hpp" // for the sparse space dense vector
#endif                                    // KRATOS_USE_AMATRIX

// External includes
#include "pybind11/pybind11.h"

// Project includes
#include "includes/model_part.h"
#include "linear_solvers/linear_solver.h"
#include "spaces/ublas_space.h"

#include "Epetra_FEVector.h"
#include "trilinos_space.h"

#include "add_trilinos_auxiliary_processes_to_python.h"

// Application includes
#include "custom_processes/auxiliary_processes/rans_wall_distance_calculation_process.h"

namespace Kratos
{

namespace Python
{
void AddTrilinosAuxiliaryProcessesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    using LocalSpaceType = UblasSpace<double, Matrix, Vector>;
    using MPISparseSpaceType = TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>;
    using MPILinearSolverType = LinearSolver<MPISparseSpaceType, LocalSpaceType>;

    using MPIRansWallDistanceCalculationProcessType = RansWallDistanceCalculationProcess<MPISparseSpaceType, LocalSpaceType, MPILinearSolverType>;
    py::class_<MPIRansWallDistanceCalculationProcessType, MPIRansWallDistanceCalculationProcessType::Pointer, Process>(
        m, "MPIRansWallDistanceCalculationProcess")
        .def(py::init<Model&, Parameters&>())
        .def("SetBuilderAndSolver", &MPIRansWallDistanceCalculationProcessType::SetBuilderAndSolver)
        .def("GetLinearSolver", &MPIRansWallDistanceCalculationProcessType::GetLinearSolver);
}

} // namespace Python.
} // Namespace Kratos
