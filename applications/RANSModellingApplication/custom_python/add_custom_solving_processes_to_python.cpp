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

#ifdef KRATOS_USING_MPI // mpi-parallel compilation
#include "Epetra_FEVector.h"
#include "trilinos_space.h"
#endif

#include "custom_python/add_custom_solving_processes_to_python.h"

// Application includes
#include "custom_processes/solving_strategies/k_epsilon_co_solving_process.h"
#include "custom_processes/solving_strategies/scalar_co_solving_process.h"

namespace Kratos
{
namespace Python
{
void AddCustomSolvingProcessesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType> LinearSolverType;

    // Adding solving strategies
    typedef ScalarCoSolvingProcess<SparseSpaceType, LocalSpaceType, LinearSolverType> ScalarCoSolvingProcessType;
    py::class_<ScalarCoSolvingProcessType, ScalarCoSolvingProcessType::Pointer, Process>(
        m, "ScalarCoSolvingProcess")
        .def(py::init<ModelPart&, Parameters&, Variable<double>&>())
        .def("AddStrategy", &ScalarCoSolvingProcessType::AddStrategy)
        .def("AddAuxiliaryProcess", &ScalarCoSolvingProcessType::AddAuxiliaryProcess);

    typedef KEpsilonCoSolvingProcess<SparseSpaceType, LocalSpaceType, LinearSolverType> KEpsilonCoSolvingProcessType;
    py::class_<KEpsilonCoSolvingProcessType, KEpsilonCoSolvingProcessType::Pointer, ScalarCoSolvingProcessType, Process>(
        m, "KEpsilonCoSolvingProcess")
        .def(py::init<ModelPart&, Parameters&>());

#ifdef KRATOS_USING_MPI // mpi-parallel compilation
    typedef TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector> MPISparseSpaceType;
    typedef LinearSolver<MPISparseSpaceType, LocalSpaceType> MPILinearSolverType;

    typedef ScalarCoSolvingProcess<MPISparseSpaceType, LocalSpaceType, MPILinearSolverType> MPIScalarCoSolvingProcessType;
    py::class_<MPIScalarCoSolvingProcessType, MPIScalarCoSolvingProcessType::Pointer, Process>(
        m, "MPIScalarCoSolvingProcess")
        .def(py::init<ModelPart&, Parameters&, Variable<double>&>())
        .def("AddStrategy", &MPIScalarCoSolvingProcessType::AddStrategy)
        .def("AddAuxiliaryProcess", &MPIScalarCoSolvingProcessType::AddAuxiliaryProcess);

    typedef KEpsilonCoSolvingProcess<MPISparseSpaceType, LocalSpaceType, MPILinearSolverType> MPIKEpsilonCoSolvingProcessType;
    py::class_<MPIKEpsilonCoSolvingProcessType, MPIKEpsilonCoSolvingProcessType::Pointer, MPIScalarCoSolvingProcessType, Process>(
        m, "MPIKEpsilonCoSolvingProcess")
        .def(py::init<ModelPart&, Parameters&>());
#endif
}

} // namespace Python.
} // Namespace Kratos
