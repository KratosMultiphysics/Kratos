//  KRATOS  _____     _ _ _
//         |_   _| __(_) (_)_ __   ___  ___
//           | || '__| | | | '_ \ / _ \/ __|
//           | || |  | | | | | | | (_) \__
//           |_||_|  |_|_|_|_| |_|\___/|___/ APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

#if defined(KRATOS_PYTHON)

// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "trilinos_space.h"
#include "spaces/ublas_space.h"
#include "includes/model_part.h"

// Linear solvers
#include "linear_solvers/linear_solver.h"
#include "linear_solvers/fallback_linear_solver.h"
#include "custom_factories/trilinos_linear_solver_factory.h"

//teuchos parameter list
#include "Teuchos_ParameterList.hpp"

#ifndef TRILINOS_EXCLUDE_AZTEC_SOLVER
#include "external_includes/aztec_solver.h"
#endif

#ifndef TRILINOS_EXCLUDE_AMESOS_SOLVER
#include "external_includes/amesos_solver.h"
#endif

#ifndef TRILINOS_EXCLUDE_AMESOS2_SOLVER
#include "external_includes/amesos2_solver.h"
#endif

#ifndef TRILINOS_EXCLUDE_ML_SOLVER
#include "external_includes/ml_solver.h"
#endif

#include "external_includes/amgcl_mpi_solver.h"
#include "external_includes/amgcl_mpi_schur_complement_solver.h"

namespace Kratos::Python
{

namespace py = pybind11;

using TrilinosSparseSpaceType = TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>;
using TrilinosLocalSpaceType = UblasSpace<double, Matrix, Vector>;
using TrilinosLinearSolverType = LinearSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType>;

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
        .def("Solve", Solve)
        ;

    using TrilinosFallbackLinearSolverType = FallbackLinearSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType>;
    py::class_<TrilinosFallbackLinearSolverType, TrilinosFallbackLinearSolverType::Pointer, TrilinosLinearSolverType>(m, "TrilinosFallbackLinearSolver")
        .def(py::init<Parameters>())
        .def(py::init<const std::vector<TrilinosLinearSolverType::Pointer>&, Parameters>())
        // .def("AddSolver", [](TrilinosFallbackLinearSolverType& rSelf, TrilinosLinearSolverType::Pointer pSolver) {
        //     rSelf.AddSolver(pSolver);
        // })
        // .def("AddSolver", [](TrilinosFallbackLinearSolverType& rSelf, const Parameters ThisParameters) {
        //     rSelf.AddSolver(ThisParameters);
        // })
        .def("GetSolvers", &TrilinosFallbackLinearSolverType::GetSolvers)
        // .def("SetSolvers", &TrilinosFallbackLinearSolverType::SetSolvers)
        .def("GetResetSolverEachTry", &TrilinosFallbackLinearSolverType::GetResetSolverEachTry)
        // .def("SetResetSolverIndexEachTry", &TrilinosFallbackLinearSolverType::SetResetSolverIndexEachTry)
        .def("GetParameters", &TrilinosFallbackLinearSolverType::GetParameters)
        .def("GetCurrentSolverIndex", &TrilinosFallbackLinearSolverType::GetCurrentSolverIndex)
        .def("ClearCurrentSolverIndex", &TrilinosFallbackLinearSolverType::ClearCurrentSolverIndex)
        ;

#ifndef TRILINOS_EXCLUDE_AZTEC_SOLVER
    typedef AztecSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType > AztecSolverType;
    py::class_<AztecSolverType, typename AztecSolverType::Pointer, TrilinosLinearSolverType >
    (m,"AztecSolver")
        .def(py::init< Teuchos::ParameterList&, std::string, Teuchos::ParameterList&, double, int, int >())
        .def(py::init<Parameters>())
        .def("SetScalingType", &AztecSolverType::SetScalingType)
        .def("__str__", PrintObject<AztecSolverType>)
        ;

    py::enum_<AztecScalingType>(m,"AztecScalingType")
        .value("NoScaling", NoScaling)
        .value("LeftScaling", LeftScaling)
        .value("SymmetricScaling", SymmetricScaling)
        ;
#endif

#ifndef TRILINOS_EXCLUDE_AMESOS_SOLVER
    typedef AmesosSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType > AmesosSolverType;
    py::class_<AmesosSolverType, typename AmesosSolverType::Pointer, TrilinosLinearSolverType >
    (m,"AmesosSolver").def( py::init<const std::string&, Teuchos::ParameterList& >())
        .def(py::init<Parameters>())
        .def_static("HasSolver", &AmesosSolverType::HasSolver)
        .def("__str__", PrintObject<AmesosSolverType>)
        ;
#endif

#ifndef TRILINOS_EXCLUDE_AMESOS2_SOLVER
    typedef Amesos2Solver<TrilinosSparseSpaceType, TrilinosLocalSpaceType > Amesos2SolverType;
    py::class_<Amesos2SolverType, typename Amesos2SolverType::Pointer, TrilinosLinearSolverType >
    (m,"Amesos2Solver").def( py::init<const std::string&, Teuchos::ParameterList& >())
        .def(py::init<Parameters>())
        .def_static("HasSolver", &Amesos2SolverType::HasSolver)
        .def("__str__", PrintObject<Amesos2SolverType>)
        ;
#endif

#ifndef TRILINOS_EXCLUDE_ML_SOLVER
    typedef MultiLevelSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType > MLSolverType;
    py::class_<MLSolverType, typename MLSolverType::Pointer, TrilinosLinearSolverType >
    (m,"MultiLevelSolver").def( py::init<Teuchos::ParameterList&, Teuchos::ParameterList&, double, int >())
        .def(py::init<Parameters>())
        .def("SetScalingType", &MLSolverType::SetScalingType)
        .def("SetReformPrecAtEachStep", &MLSolverType::SetReformPrecAtEachStep)
        .def("__str__", PrintObject<MLSolverType>)
        .def_static("SetDefaults", &MLSolverType::SetDefaults)
        ;

    py::enum_<MLSolverType::ScalingType>(m,"MLSolverScalingType")
        .value("NoScaling", MLSolverType::NoScaling)
        .value("LeftScaling", MLSolverType::LeftScaling)
        ;
#endif

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

    typedef LinearSolverFactory< TrilinosSparseSpaceType, TrilinosLocalSpaceType > TrilinosLinearSolverFactoryType;

    py::class_<TrilinosLinearSolverFactoryType, TrilinosLinearSolverFactoryType::Pointer>(m, "TrilinosLinearSolverFactory")
        .def( py::init< >() )
        .def("Create",&TrilinosLinearSolverFactoryType::Create)
        .def("Has",&TrilinosLinearSolverFactoryType::Has)
        ;
}

} // namespace Python:: Kratos.

#endif // KRATOS_PYTHON defined
