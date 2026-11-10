//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Raul Bravo
//

#ifdef KRATOS_USE_FUTURE

// System includes

// External includes
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

// Project includes
#include "includes/define.h"

// Application includes
#include "custom_python/add_future_utilities_to_python.h"
#include "future/rom_projector.h"
#include "future/hrom_scheme.h"
#include "future/solving_strategies/schemes/static_scheme.h"

namespace Kratos::Python {

void AddFutureUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    using RomProjectorType = Future::RomProjector<Future::SerialLinearAlgebraTraits>;
    py::class_<RomProjectorType, typename RomProjectorType::Pointer>(m, "RomProjector")
        .def(py::init<typename RomProjectorType::SchemeType::Pointer, typename RomProjectorType::StrategyDataType::Pointer>(),
            py::arg("scheme"), py::arg("strategy_data"),
            "Creates the projector from a Future implicit scheme and the strategy data it has initialized.")
        .def("BuildEffectiveSystem", &RomProjectorType::BuildEffectiveSystem,
            "Builds the full order system from the current database and applies the constraints.")
        .def("Project", [](RomProjectorType& rSelf, const Eigen::Ref<const RomProjectorType::EigenDynamicMatrix>& rPhi){
                rSelf.Project(rPhi);
                return py::make_tuple(rSelf.GetReducedLhs(), rSelf.GetReducedRhs());
            }, py::arg("phi"),
            "Galerkin projection of the effective system onto the basis phi. Returns the reduced LHS and RHS.")
        .def("SolveReduced", &RomProjectorType::SolveReduced,
            "Solves the reduced system of the last projection and returns the increment of the reduced coordinates.")
        .def("SetSolution", &RomProjectorType::SetSolution, py::arg("solution"),
            "Sets the provided solution vector in the free DOFs.")
        ;

    using StaticSchemeType = Future::StaticScheme<Future::SerialLinearAlgebraTraits>;
    using HRomStaticSchemeType = Future::HRomScheme<StaticSchemeType>;
    py::class_<HRomStaticSchemeType, typename HRomStaticSchemeType::Pointer, StaticSchemeType>(m, "HRomStaticScheme")
        .def(py::init<ModelPart&, Parameters>(), py::arg("model_part"), py::arg("settings"),
            "Static scheme that only assembles the HROM entities, scaled by their weights.")
        .def("SetHRomEntities", &HRomStaticSchemeType::SetHRomEntities,
            py::arg("element_ids"), py::arg("element_weights"), py::arg("condition_ids"), py::arg("condition_weights"),
            "Sets the HROM elements and conditions and their weights (the ones with a null weight are skipped).")
        .def("NumberOfHRomElements", &HRomStaticSchemeType::NumberOfHRomElements)
        .def("NumberOfHRomConditions", &HRomStaticSchemeType::NumberOfHRomConditions)
        ;
}

}  // namespace Kratos::Python.

#endif // KRATOS_USE_FUTURE
