//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

// System includes


// External includes


// Project includes
#include "add_custom_utilities_to_python.h"
#include "custom_utilities/estimate_dt_utility.h"
#include "custom_utilities/shallow_water_utilities.h"
#include "custom_utilities/move_shallow_mesh_utility.h"
#include "custom_utilities/derivatives_recovery_utility.h"
#include "custom_utilities/interpolate_sw_to_pfem_utility.hpp"


namespace Kratos
{

namespace Python
{

typedef ModelPart::NodesContainerType NodesContainerType;

typedef ModelPart::ElementsContainerType ElementsContainerType;

typedef ModelPart::ConditionsContainerType ConditionsContainerType;

typedef ModelPart::PropertiesContainerType PropertiesContainerType;

template<class TContainerType>
array_1d<double,3> ComputeHydrostaticForces1(
    ShallowWaterUtilities& rUtility,
    TContainerType& rContainer,
    const ProcessInfo& rProcessInfo)
{
    return rUtility.ComputeHydrostaticForces(rContainer, rProcessInfo);
}

template<class TContainerType>
array_1d<double,3> ComputeHydrostaticForces2(
    ShallowWaterUtilities& rUtility,
    TContainerType& rContainer,
    const ProcessInfo& rProcessInfo,
    const double RelativeDryHeight)
{
    return rUtility.ComputeHydrostaticForces(rContainer, rProcessInfo, RelativeDryHeight);
}

void  AddCustomUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_< ShallowWaterUtilities > (m, "ShallowWaterUtilities")
        .def(py::init<>())
        .def("ComputeFreeSurfaceElevation", &ShallowWaterUtilities::ComputeFreeSurfaceElevation)
        .def("ComputeHeightFromFreeSurface", &ShallowWaterUtilities::ComputeHeightFromFreeSurface)
        .def("ComputeVelocity", &ShallowWaterUtilities::ComputeVelocity)
        .def("ComputeMomentum", &ShallowWaterUtilities::ComputeMomentum)
        .def("ComputeFroude", &ShallowWaterUtilities::ComputeFroude<true>)
        .def("ComputeFroudeNonHistorical", &ShallowWaterUtilities::ComputeFroude<false>)
        .def("ComputeEnergy", &ShallowWaterUtilities::ComputeEnergy<true>)
        .def("ComputeEnergyNonHistorical", &ShallowWaterUtilities::ComputeEnergy<false>)
        .def("FlipScalarVariable", &ShallowWaterUtilities::FlipScalarVariable)
        .def("IdentifySolidBoundary", &ShallowWaterUtilities::IdentifySolidBoundary)
        .def("NormalizeVector", &ShallowWaterUtilities::NormalizeVector)
        .def("SmoothHistoricalVariable", &ShallowWaterUtilities::SmoothHistoricalVariable<double>)
        .def("SmoothHistoricalVariable", &ShallowWaterUtilities::SmoothHistoricalVariable<array_1d<double,3>>)
        .def("CopyVariableToPreviousTimeStep", &ShallowWaterUtilities::CopyVariableToPreviousTimeStep<Variable<double>&>)
        .def("CopyVariableToPreviousTimeStep", &ShallowWaterUtilities::CopyVariableToPreviousTimeStep<Variable<array_1d<double,3>>&>)
        .def("SetMinimumValue", &ShallowWaterUtilities::SetMinimumValue)
        .def("SetMeshZCoordinateToZero", &ShallowWaterUtilities::SetMeshZCoordinateToZero)
        .def("SetMeshZ0CoordinateToZero", &ShallowWaterUtilities::SetMeshZ0CoordinateToZero)
        .def("SetMeshZCoordinate", &ShallowWaterUtilities::SetMeshZCoordinate)
        .def("OffsetMeshZCoordinate", &ShallowWaterUtilities::OffsetMeshZCoordinate)
        .def("SwapYZCoordinates", &ShallowWaterUtilities::SwapYZCoordinates)
        .def("SwapY0Z0Coordinates", &ShallowWaterUtilities::SwapY0Z0Coordinates)
        .def("SwapYZComponents", &ShallowWaterUtilities::SwapYZComponents)
        .def("SwapYZComponentsNonHistorical", &ShallowWaterUtilities::SwapYZComponentsNonHistorical<NodesContainerType>)
        .def("SwapYZComponentsNonHistorical", &ShallowWaterUtilities::SwapYZComponentsNonHistorical<ElementsContainerType>)
        .def("SwapYZComponentsNonHistorical", &ShallowWaterUtilities::SwapYZComponentsNonHistorical<ConditionsContainerType>)
        .def("StoreNonHistoricalGiDNoDataIfDry", &ShallowWaterUtilities::StoreNonHistoricalGiDNoDataIfDry)
        .def("ComputeL2Norm", &ShallowWaterUtilities::ComputeL2Norm<true>)
        .def("ComputeL2Norm", &ShallowWaterUtilities::ComputeL2NormAABB<true>)
        .def("ComputeL2NormNonHistorical", &ShallowWaterUtilities::ComputeL2Norm<false>)
        .def("ComputeL2NormNonHistorical", &ShallowWaterUtilities::ComputeL2NormAABB<false>)
        .def("ComputeHydrostaticForces", ComputeHydrostaticForces1<ElementsContainerType>)
        .def("ComputeHydrostaticForces", ComputeHydrostaticForces2<ElementsContainerType>)
        .def("ComputeHydrostaticForces", ComputeHydrostaticForces1<ConditionsContainerType>)
        .def("ComputeHydrostaticForces", ComputeHydrostaticForces2<ConditionsContainerType>)
        .def("OffsetIds", [](ShallowWaterUtilities& self, NodesContainerType&      rContainer){self.OffsetIds(rContainer);})
        .def("OffsetIds", [](ShallowWaterUtilities& self, ElementsContainerType&   rContainer){self.OffsetIds(rContainer);})
        .def("OffsetIds", [](ShallowWaterUtilities& self, ConditionsContainerType& rContainer){self.OffsetIds(rContainer);})
        .def("OffsetIds", [](ShallowWaterUtilities& self, PropertiesContainerType& rContainer){self.OffsetIds(rContainer);})
        .def("OffsetIds", [](ShallowWaterUtilities& self, NodesContainerType&      rContainer, const double Value){self.OffsetIds(rContainer, Value);})
        .def("OffsetIds", [](ShallowWaterUtilities& self, ElementsContainerType&   rContainer, const double Value){self.OffsetIds(rContainer, Value);})
        .def("OffsetIds", [](ShallowWaterUtilities& self, ConditionsContainerType& rContainer, const double Value){self.OffsetIds(rContainer, Value);})
        .def("OffsetIds", [](ShallowWaterUtilities& self, PropertiesContainerType& rContainer, const double Value){self.OffsetIds(rContainer, Value);})
        ;

    py::class_< EstimateTimeStepUtility > (m, "EstimateTimeStepUtility")
        .def(py::init<ModelPart&, Parameters>())
        .def("Execute", &EstimateTimeStepUtility::Execute)
        ;

    py::class_<MoveShallowMeshUtility>(m, "MoveShallowMeshUtility")
        .def(py::init<ModelPart&, ModelPart&, Parameters>())
        .def("Check", &MoveShallowMeshUtility::Check)
        .def("Initialize", &MoveShallowMeshUtility::Initialize)
        .def("MoveMesh", &MoveShallowMeshUtility::MoveMesh)
        .def("MapResults", &MoveShallowMeshUtility::MapResults)
        ;

    py::class_<DerivativesRecoveryUtility<2>>(m, "DerivativesRecoveryUtility2D")
        .def_static("Check", &DerivativesRecoveryUtility<2>::Check)
        .def_static("CalculatePolynomialWeights", &DerivativesRecoveryUtility<2>::CalculatePolynomialWeights)
        .def_static("RecoverDivergence", &DerivativesRecoveryUtility<2>::RecoverDivergence)
        .def_static("RecoverGradient", &DerivativesRecoveryUtility<2>::RecoverGradient)
        .def_static("RecoverLaplacian", [](ModelPart& rModelPart,
            const Variable<double>& rOriginVariable,
            const Variable<double>& rDestinationVariable,
            const std::size_t BufferStep) {
                DerivativesRecoveryUtility<2>::RecoverLaplacian(rModelPart, rOriginVariable, rDestinationVariable, BufferStep);})
        .def_static("RecoverLaplacian", [](ModelPart& rModelPart,
            const Variable<array_1d<double,3>>& rOriginVariable,
            const Variable<array_1d<double,3>>& rDestinationVariable,
            const std::size_t BufferStep) {
                DerivativesRecoveryUtility<2>::RecoverLaplacian(rModelPart, rOriginVariable, rDestinationVariable, BufferStep);})
        ;

    py::class_<DerivativesRecoveryUtility<3>>(m, "DerivativesRecoveryUtility3D")
        .def_static("Check", &DerivativesRecoveryUtility<3>::Check)
        .def_static("CalculatePolynomialWeights", &DerivativesRecoveryUtility<3>::CalculatePolynomialWeights)
        .def_static("RecoverDivergence", &DerivativesRecoveryUtility<3>::RecoverDivergence)
        .def_static("RecoverGradient", &DerivativesRecoveryUtility<3>::RecoverGradient)
        .def_static("RecoverLaplacian", [](ModelPart& rModelPart,
            const Variable<double>& rOriginVariable,
            const Variable<double>& rDestinationVariable,
            const std::size_t BufferStep) {
                DerivativesRecoveryUtility<3>::RecoverLaplacian(rModelPart, rOriginVariable, rDestinationVariable, BufferStep);})
        .def_static("RecoverLaplacian", [](ModelPart& rModelPart,
            const Variable<array_1d<double,3>>& rOriginVariable,
            const Variable<array_1d<double,3>>& rDestinationVariable,
            const std::size_t BufferStep) {
                DerivativesRecoveryUtility<3>::RecoverLaplacian(rModelPart, rOriginVariable, rDestinationVariable, BufferStep);})
        ;

    py::class_< InterpolateSwToPfemUtility >
        (m, "InterpolateSwToPfemUtility")
        .def( py::init<>())
        .def("InterpolateVariables",&InterpolateSwToPfemUtility::InterpolateVariables)
    ;

}

}  // namespace Python.

} // Namespace Kratos
