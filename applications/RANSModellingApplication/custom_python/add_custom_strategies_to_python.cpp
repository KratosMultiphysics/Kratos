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
void AddCustomStrategiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    using LocalSpaceType = UblasSpace<double, Matrix, Vector>;
    using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
    using BaseSchemeType = Scheme<SparseSpaceType, LocalSpaceType>;

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
    using MPISparseSpaceType = TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>;

    using MPIBaseSchemeType = Scheme<MPISparseSpaceType, LocalSpaceType>;

    using MPIConvergenceCriteria = ConvergenceCriteria<MPISparseSpaceType, LocalSpaceType>;

    py::class_<MPIGenericConvergenceCriteria<MPISparseSpaceType, LocalSpaceType>,
               typename MPIGenericConvergenceCriteria<MPISparseSpaceType, LocalSpaceType>::Pointer, MPIConvergenceCriteria>(
        m, "MPIGenericScalarConvergenceCriteria")
        .def(py::init<MPISparseSpaceType::DataType, MPISparseSpaceType::DataType>());

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
