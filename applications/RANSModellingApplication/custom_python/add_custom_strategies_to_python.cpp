//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

// System includes

// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "custom_python/add_custom_strategies_to_python.h"
#include "includes/define_python.h"
#include "spaces/ublas_space.h"

#ifdef KRATOS_USING_MPI // mpi-parallel compilation
#include "Epetra_FEVector.h"
#include "trilinos_space.h"
#endif

// strategies
#include "custom_strategies/generic_residual_based_bossak_velocity_scalar_scheme.h"
#include "custom_strategies/generic_residualbased_simple_steady_scalar_scheme.h"

// convergence criterians
#include "custom_strategies/generic_convergence_criteria.h"
#ifdef KRATOS_USING_MPI
#include "custom_strategies/generic_convergence_criteria_mpi.h"
#endif

namespace Kratos
{
namespace Python
{
typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
#ifdef KRATOS_USING_MPI
typedef TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector> MPISparseSpaceType;
void MoveMesh(Scheme<MPISparseSpaceType, LocalSpaceType>& dummy,
              ModelPart::NodesContainerType& rNodes)
{
    for (ModelPart::NodeIterator i = rNodes.begin(); i != rNodes.end(); ++i)
    {
        const array_1d<double, 3>& disp = i->FastGetSolutionStepValue(DISPLACEMENT);
        (i)->X() = (i)->X0() + disp[0];
        (i)->Y() = (i)->Y0() + disp[1];
        (i)->Z() = (i)->Z0() + disp[2];
    }
}
#endif

void AddCustomStrategiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef Scheme<SparseSpaceType, LocalSpaceType> BaseSchemeType;

    // Convergence criteria
    py::class_<GenericConvergenceCriteria<SparseSpaceType, LocalSpaceType>,
               typename GenericConvergenceCriteria<SparseSpaceType, LocalSpaceType>::Pointer,
               ConvergenceCriteria<SparseSpaceType, LocalSpaceType>>(
        m, "GenericScalarConvergenceCriteria")
        .def(py::init<double, double>());

    py::class_<GenericResidualBasedBossakVelocityScalarScheme<SparseSpaceType, LocalSpaceType>,
               typename GenericResidualBasedBossakVelocityScalarScheme<SparseSpaceType, LocalSpaceType>::Pointer, BaseSchemeType>(
        m, "GenericResidualBasedBossakVelocityDynamicScalarScheme")
        .def(py::init<const double, const double, const Variable<double>&,
                      const Variable<double>&, const Variable<double>&>());

    py::class_<GenericResidualBasedSimpleSteadyScalarScheme<SparseSpaceType, LocalSpaceType>,
               typename GenericResidualBasedSimpleSteadyScalarScheme<SparseSpaceType, LocalSpaceType>::Pointer, BaseSchemeType>(
        m, "GenericResidualBasedSimpleSteadyScalarScheme")
        .def(py::init<const double>());

#ifdef KRATOS_USING_MPI // mpi-parallel compilation

    typedef Scheme<MPISparseSpaceType, LocalSpaceType> MPIBaseSchemeType;

    typedef ConvergenceCriteria<MPISparseSpaceType, LocalSpaceType> MPIConvergenceCriteria;

    py::class_<MPIConvergenceCriteria, MPIConvergenceCriteria::Pointer>(
        m, "MPIConvergenceCriteria")
        .def(py::init<>())
        .def("SetActualizeRHSFlag", &MPIConvergenceCriteria::SetActualizeRHSFlag)
        .def("GetActualizeRHSflag", &MPIConvergenceCriteria::GetActualizeRHSflag)
        .def("PreCriteria", &MPIConvergenceCriteria::PreCriteria)
        .def("PostCriteria", &MPIConvergenceCriteria::PostCriteria)
        .def("Initialize", &MPIConvergenceCriteria::Initialize)
        .def("InitializeSolutionStep", &MPIConvergenceCriteria::InitializeSolutionStep)
        .def("FinalizeSolutionStep", &MPIConvergenceCriteria::FinalizeSolutionStep)
        .def("Check", &MPIConvergenceCriteria::Check)
        .def("SetEchoLevel", &MPIConvergenceCriteria::SetEchoLevel);

    py::class_<MPIGenericConvergenceCriteria<MPISparseSpaceType, LocalSpaceType>,
               typename MPIGenericConvergenceCriteria<MPISparseSpaceType, LocalSpaceType>::Pointer, MPIConvergenceCriteria>(
        m, "MPIGenericScalarConvergenceCriteria")
        .def(py::init<MPISparseSpaceType::DataType, MPISparseSpaceType::DataType>());

    py::class_<MPIBaseSchemeType, typename MPIBaseSchemeType::Pointer>(
        m, "MPIScheme")
        .def(py::init<>())
        .def("Initialize", &MPIBaseSchemeType::Initialize)
        .def("SchemeIsInitialized", &MPIBaseSchemeType::SchemeIsInitialized)
        .def("ElementsAreInitialized", &MPIBaseSchemeType::ElementsAreInitialized)
        .def("ConditionsAreInitialized", &MPIBaseSchemeType::ConditionsAreInitialized)
        .def("InitializeElements", &MPIBaseSchemeType::InitializeElements)
        .def("InitializeConditions", &MPIBaseSchemeType::InitializeConditions)
        .def("InitializeSolutionStep", &MPIBaseSchemeType::InitializeSolutionStep)
        .def("FinalizeSolutionStep", &MPIBaseSchemeType::FinalizeSolutionStep)
        .def("InitializeNonLinIteration", &MPIBaseSchemeType::InitializeNonLinIteration)
        .def("FinalizeNonLinIteration", &MPIBaseSchemeType::FinalizeNonLinIteration)
        .def("Predict", &MPIBaseSchemeType::Predict)
        .def("Update", &MPIBaseSchemeType::Update)
        .def("CalculateOutputData", &MPIBaseSchemeType::CalculateOutputData)
        .def("Clean", &MPIBaseSchemeType::Clean)
        .def("Clear", &MPIBaseSchemeType::Clear)
        .def("MoveMesh", MoveMesh)
        .def("Check", &MPIBaseSchemeType::Check);

    py::class_<GenericResidualBasedBossakVelocityScalarScheme<MPISparseSpaceType, LocalSpaceType>,
               typename GenericResidualBasedBossakVelocityScalarScheme<MPISparseSpaceType, LocalSpaceType>::Pointer, MPIBaseSchemeType>(
        m, "MPIGenericResidualBasedBossakVelocityDynamicScalarScheme")
        .def(py::init<const double, const double, const Variable<double>&,
                      const Variable<double>&, const Variable<double>&>());

    py::class_<GenericResidualBasedSimpleSteadyScalarScheme<MPISparseSpaceType, LocalSpaceType>,
               typename GenericResidualBasedSimpleSteadyScalarScheme<MPISparseSpaceType, LocalSpaceType>::Pointer, MPIBaseSchemeType>(
        m, "MPIGenericResidualBasedSimpleSteadyScalarScheme")
        .def(py::init<const double>());
#endif
}

} // namespace Python.
} // Namespace Kratos
