//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//


// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/hydraulic_fluid_auxiliary_utilities.h"

namespace Kratos::Python {

void AddCustomUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;
    py::class_<HydraulicFluidAuxiliaryUtilities>(m, "HydraulicFluidAuxiliaryUtilities")
        .def_static("CalculateWettedPetimeter", [](ModelPart &rModelPart, const Flags &rSkinFlag, const Variable<double> &rDistanceVariable, const bool IsHistorical)
                    { return HydraulicFluidAuxiliaryUtilities::CalculateWettedPetimeter(rModelPart, rSkinFlag, rDistanceVariable, IsHistorical); })
        .def_static("CalculateWettedArea", [](ModelPart &rModelPart, const Flags &rSkinFlag, const Variable<double> &rDistanceVariable, const bool IsHistorical)
                    { return HydraulicFluidAuxiliaryUtilities::CalculateWettedArea(rModelPart, rSkinFlag, rDistanceVariable, IsHistorical); })
        .def_static("InitialWaterDepth", [](ModelPart &rModelPart)
                    { return HydraulicFluidAuxiliaryUtilities::InitialWaterDepth(rModelPart); })
        .def_static("SetInletVelocity", [](ModelPart &rModelPart, double InletVelocity, const Variable<double> &rDistanceVariable)
                    { return HydraulicFluidAuxiliaryUtilities::SetInletVelocity(rModelPart, InletVelocity, rDistanceVariable); })
        .def_static("FreeInlet", [](ModelPart &rModelPart)
                    { return HydraulicFluidAuxiliaryUtilities::FreeInlet(rModelPart); })
        .def_static("SetInletFreeSurface", [](ModelPart &rModelPart, const Flags &rSkinFlag, const Variable<double> &rDistanceVariable)
                    { return HydraulicFluidAuxiliaryUtilities::SetInletFreeSurface(rModelPart, rSkinFlag, rDistanceVariable); })
        .def_static("FixCornerNodeVelocity", [](ModelPart &rModelPart, double MaximumAngle)
                    { return HydraulicFluidAuxiliaryUtilities::FixCornerNodeVelocity(rModelPart, MaximumAngle); })
        .def_static("TurnOffGravityOnAirElements", [](ModelPart &rModelPart)
                    { return HydraulicFluidAuxiliaryUtilities::TurnOffGravityOnAirElements(rModelPart); })
        .def_static("MaximumWaterDepthChange", [](ModelPart &rModelPart)
                    { return HydraulicFluidAuxiliaryUtilities::MaximumWaterDepthChange(rModelPart); })
        .def_static("CalculateArtificialViscosity", [](ModelPart &rModelPart, double LimiterCoefficient)
                    { return HydraulicFluidAuxiliaryUtilities::CalculateArtificialViscosity(rModelPart, LimiterCoefficient); });
}

} // namespace Kratos::Python
