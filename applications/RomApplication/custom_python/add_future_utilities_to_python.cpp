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
#include "future/legacy_scheme_adapter.h"
#include "future/entity_residuals_scheme.h"
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
        .def("ProjectEntityResiduals", &RomProjectorType::ProjectEntityResiduals, py::arg("phi"),
            "Projects the residual of each element and condition (elements first) onto the basis phi, one row per entity.")
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

    // Prototype of the ROM as the solution strategy of the solvers (see python_scripts/future/future_rom_strategy.py)
    using LegacySchemeAdapterType = Future::LegacySchemeAdapter<Future::SerialLinearAlgebraTraits>;
    py::class_<LegacySchemeAdapterType, typename LegacySchemeAdapterType::Pointer, StaticSchemeType>(m, "LegacySchemeAdapter")
        .def(py::init<ModelPart&, Parameters, typename LegacySchemeAdapterType::LegacySchemeType::Pointer>(),
            py::arg("model_part"), py::arg("settings"), py::arg("legacy_scheme"),
            "Future scheme (DOFs, system and assembly) whose time integration is the one of the given legacy scheme.")
        ;

    // Schemes of the ROM (not HROM) runs: the base scheme plus the projection of the entity residuals of the HROM training
    using StrategyDataType = Future::ImplicitStrategyData<Future::SerialLinearAlgebraTraits>;
    using RomStaticSchemeType = Future::EntityResidualsScheme<StaticSchemeType>;
    py::class_<RomStaticSchemeType, typename RomStaticSchemeType::Pointer, StaticSchemeType>(m, "RomStaticScheme")
        .def(py::init<ModelPart&, Parameters>(), py::arg("model_part"), py::arg("settings"))
        .def("ProjectEntityResiduals", [](RomStaticSchemeType& rSelf, const Eigen::Ref<const RomStaticSchemeType::EigenDynamicMatrix>& rPhi, StrategyDataType& rData){
            return rSelf.ProjectEntityResiduals(rPhi, rData);}, py::arg("phi"), py::arg("strategy_data"))
        ;
    using RomLegacySchemeAdapterType = Future::EntityResidualsScheme<LegacySchemeAdapterType>;
    py::class_<RomLegacySchemeAdapterType, typename RomLegacySchemeAdapterType::Pointer, LegacySchemeAdapterType>(m, "RomLegacySchemeAdapter")
        .def(py::init<ModelPart&, Parameters, typename LegacySchemeAdapterType::LegacySchemeType::Pointer>(),
            py::arg("model_part"), py::arg("settings"), py::arg("legacy_scheme"))
        .def("ProjectEntityResiduals", [](RomLegacySchemeAdapterType& rSelf, const Eigen::Ref<const RomLegacySchemeAdapterType::EigenDynamicMatrix>& rPhi, StrategyDataType& rData){
            return rSelf.ProjectEntityResiduals(rPhi, rData);}, py::arg("phi"), py::arg("strategy_data"))
        ;

    using HRomLegacySchemeAdapterType = Future::HRomScheme<LegacySchemeAdapterType>;
    py::class_<HRomLegacySchemeAdapterType, typename HRomLegacySchemeAdapterType::Pointer, LegacySchemeAdapterType>(m, "HRomLegacySchemeAdapter")
        .def(py::init<ModelPart&, Parameters, typename LegacySchemeAdapterType::LegacySchemeType::Pointer>(),
            py::arg("model_part"), py::arg("settings"), py::arg("legacy_scheme"),
            "LegacySchemeAdapter that only assembles the HROM entities, scaled by their weights.")
        .def("SetHRomEntities", &HRomLegacySchemeAdapterType::SetHRomEntities,
            py::arg("element_ids"), py::arg("element_weights"), py::arg("condition_ids"), py::arg("condition_weights"))
        .def("NumberOfHRomElements", &HRomLegacySchemeAdapterType::NumberOfHRomElements)
        .def("NumberOfHRomConditions", &HRomLegacySchemeAdapterType::NumberOfHRomConditions)
        ;
}

}  // namespace Kratos::Python.

#endif // KRATOS_USE_FUTURE
