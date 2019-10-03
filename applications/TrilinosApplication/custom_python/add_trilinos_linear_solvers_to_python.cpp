//  KRATOS  _____     _ _ _
//         |_   _| __(_) (_)_ __   ___  ___
//           | || '__| | | | '_ \ / _ \/ __|
//           | || |  | | | | | | | (_) \__
//           |_||_|  |_|_|_|_| |_|\___/|___/ APPLICATION
//
//  License:             BSD License
//                                       Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes

#if defined(KRATOS_PYTHON)
// External includes

// Trilinos includes
#include "Epetra_FEVector.h"

// Project includes
#include "includes/define_python.h"
#include "trilinos_space.h"
#include "spaces/ublas_space.h"
#include "includes/model_part.h"

// Linear solvers
#include "linear_solvers/linear_solver.h"
#include "custom_factories/trilinos_linear_solver_factory.h"

//teuchos parameter list
#include "Teuchos_ParameterList.hpp"

#include "external_includes/aztec_solver.h"
#include "external_includes/amesos_solver.h"
#include "external_includes/ml_solver.h"

#include "external_includes/amgcl_mpi_solver.h"
#include "external_includes/amgcl_mpi_schur_complement_solver.h"

namespace Kratos
{

namespace Python
{

namespace py = pybind11;

typedef TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector> TrilinosSparseSpaceType;
typedef UblasSpace<double, Matrix, Vector> TrilinosLocalSpaceType;
typedef LinearSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType > TrilinosLinearSolverType;

void Solve(TrilinosLinearSolverType& solver,
           TrilinosSparseSpaceType::MatrixType& rA,
           TrilinosSparseSpaceType::VectorType& rX,
           TrilinosSparseSpaceType::VectorType& rB)
{
    solver.Solve(rA,rX,rB);
}

void  AddLinearSolvers(pybind11::module& m)
{
    py::class_<TrilinosLinearSolverType, TrilinosLinearSolverType::Pointer > (m,"TrilinosLinearSolver")
        .def(py::init<>())
        .def("Solve", Solve);

    typedef AztecSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType > AztecSolverType;
    py::class_<AztecSolverType, typename AztecSolverType::Pointer, TrilinosLinearSolverType >
    (m,"AztecSolver")
        .def(py::init< Teuchos::ParameterList&, std::string, Teuchos::ParameterList&, double, int, int >())
        .def(py::init<Parameters>())
        .def("SetScalingType", &AztecSolverType::SetScalingType)
        .def("__str__", PrintObject<AztecSolverType>)
        ;

    typedef AmesosSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType > AmesosSolverType;
    py::class_<AmesosSolverType, typename AmesosSolverType::Pointer, TrilinosLinearSolverType >
    (m,"AmesosSolver").def( py::init<const std::string&, Teuchos::ParameterList& >())
        .def(py::init<Parameters>())
        .def(py::init<Parameters>())
        .def_static("HasSolver", &AmesosSolverType::HasSolver)
        .def("__str__", PrintObject<AmesosSolverType>)
        ;

    typedef MultiLevelSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType > MLSolverType;
    py::class_<MLSolverType, typename MLSolverType::Pointer, TrilinosLinearSolverType >
    (m,"MultiLevelSolver").def( py::init<Teuchos::ParameterList&, Teuchos::ParameterList&, double, int >())
        .def(py::init<Parameters>())
        .def("SetScalingType", &MLSolverType::SetScalingType)
        .def("SetReformPrecAtEachStep", &MLSolverType::SetReformPrecAtEachStep)
        .def("__str__", PrintObject<MLSolverType>)
        .def_static("SetDefaults", &MLSolverType::SetDefaults)
        ;

    typedef AmgclMPISolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType > AmgclMPISolverType;
    py::class_<AmgclMPISolverType, typename AmgclMPISolverType::Pointer, TrilinosLinearSolverType >
    (m,"AmgclMPISolver")
        .def( py::init<Parameters>())
        .def("__str__", PrintObject<AmgclMPISolverType>)
        ;

    typedef AmgclMPISchurComplementSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType > AmgclMPISchurComplementSolverType;
    py::class_<AmgclMPISchurComplementSolverType, typename AmgclMPISchurComplementSolverType::Pointer, TrilinosLinearSolverType >
    (m,"AmgclMPISchurComplementSolver")
        .def( py::init<Parameters>())
        .def("__str__", PrintObject<AmgclMPISchurComplementSolverType>)
        ;

    py::enum_<AztecScalingType>(m,"AztecScalingType")
        .value("NoScaling", NoScaling)
        .value("LeftScaling", LeftScaling)
        .value("SymmetricScaling", SymmetricScaling)
        ;

    py::enum_<MLSolverType::ScalingType>(m,"MLSolverScalingType")
        .value("NoScaling", MLSolverType::NoScaling)
        .value("LeftScaling", MLSolverType::LeftScaling)
        ;

    typedef LinearSolverFactory< TrilinosSparseSpaceType, TrilinosLocalSpaceType > TrilinosLinearSolverFactoryType;

    py::class_<TrilinosLinearSolverFactoryType, TrilinosLinearSolverFactoryType::Pointer>(m, "TrilinosLinearSolverFactory")
        .def( py::init< >() )
        .def("Create",&TrilinosLinearSolverFactoryType::Create)
        .def("Has",&TrilinosLinearSolverFactoryType::Has)
        ;
}

} // namespace Python.

} // namespace Kratos.

#endif // KRATOS_PYTHON defined
