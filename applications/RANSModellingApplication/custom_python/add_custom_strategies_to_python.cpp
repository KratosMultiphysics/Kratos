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

// strategies
#include "custom_strategies/generic_residual_based_bossak_velocity_scalar_scheme.h"

// convergence criterians
#include "custom_strategies/generic_convergence_criteria.h"

namespace Kratos
{
namespace Python
{
void AddCustomStrategiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef Scheme<SparseSpaceType, LocalSpaceType> BaseSchemeType;

    // Convergence criteria
    py::class_<GenericConvergenceCriteria<SparseSpaceType, LocalSpaceType>,
               typename GenericConvergenceCriteria<SparseSpaceType, LocalSpaceType>::Pointer,
               ConvergenceCriteria<SparseSpaceType, LocalSpaceType>>(
        m, "GenericScalarConvergenceCriteria")
        .def(py::init<double, double>())
        ;

    py::class_<GenericResidualBasedBossakVelocityScalarScheme<SparseSpaceType, LocalSpaceType>,
               typename GenericResidualBasedBossakVelocityScalarScheme<SparseSpaceType, LocalSpaceType>::Pointer, BaseSchemeType>(
        m, "GenericResidualBasedBossakVelocityDynamicScalarScheme")
        .def(py::init<const double, const Variable<double>&, const Variable<double>&,
                      const Variable<double>&>())
        ;
}

} // namespace Python.
} // Namespace Kratos
