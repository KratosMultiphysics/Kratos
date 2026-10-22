//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes

// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define_python.h"
#include "spaces/default_spaces.h"

// schemes
#include "custom_strategies/bossak_relaxation_scalar_scheme.h"
#include "custom_strategies/steady_scalar_scheme.h"
#include "custom_strategies/algebraic_flux_corrected_steady_scalar_scheme.h"

// Include base h
#include "custom_python/add_custom_strategies_to_python.h"

namespace Kratos
{
namespace Python
{
void AddCustomStrategiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    using LocalSpaceType = TDefaultDenseSpace<double>;
    using SparseSpaceType = TDefaultSparseSpace<double>;
    using BaseSchemeType = Scheme<SparseSpaceType, LocalSpaceType>;

    // add schemes
    using SteadyScalarSchemeType = SteadyScalarScheme<SparseSpaceType, LocalSpaceType>;
    py::class_<SteadyScalarSchemeType, typename SteadyScalarSchemeType::Pointer, BaseSchemeType>(m, "SteadyScalarScheme")
        .def(py::init<const double>());

    using AlgebraicFluxCorrectedSteadyScalarSchemeType = AlgebraicFluxCorrectedSteadyScalarScheme<SparseSpaceType, LocalSpaceType>;
    py::class_<AlgebraicFluxCorrectedSteadyScalarSchemeType, typename AlgebraicFluxCorrectedSteadyScalarSchemeType::Pointer, BaseSchemeType>(m, "AlgebraicFluxCorrectedSteadyScalarScheme")
        .def(py::init<const double, const Flags&>())
        .def(py::init<const double, const Flags&, const Variable<int>&>());

    using BossakRelaxationScalarSchemeType = BossakRelaxationScalarScheme<SparseSpaceType, LocalSpaceType>;
    py::class_<BossakRelaxationScalarSchemeType, typename BossakRelaxationScalarSchemeType::Pointer, BaseSchemeType>(m, "BossakRelaxationScalarScheme")
        .def(py::init<const double, const double, const Variable<double>&>());

}

} // namespace Python.
} // Namespace Kratos
