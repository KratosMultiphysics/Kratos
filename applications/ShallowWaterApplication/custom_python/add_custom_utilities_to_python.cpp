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
#include "custom_utilities/move_shallow_water_particle_utility.h"
#include "custom_utilities/estimate_dt_utility.h"
#include "custom_utilities/replicate_model_part_utility.h"
#include "custom_utilities/shallow_water_utilities.h"
#include "custom_utilities/post_process_utilities.h"
#include "custom_utilities/bfecc_convection_utility.h"


namespace Kratos
{

namespace Python
{

  void  AddCustomUtilitiesToPython(pybind11::module& m)
  {
    namespace py = pybind11;

    py::class_< MoveShallowWaterParticleUtility<2> > (m, "MoveShallowWaterParticleUtility")
        .def(py::init<ModelPart& , Parameters >())
        .def("MountBin", &MoveShallowWaterParticleUtility<2>::MountBin)
        .def("MoveParticles", &MoveShallowWaterParticleUtility<2>::MoveParticles)
        .def("CorrectParticlesWithoutMovingUsingDeltaVariables", &MoveShallowWaterParticleUtility<2>::CorrectParticlesWithoutMovingUsingDeltaVariables)
        .def("PreReseed", &MoveShallowWaterParticleUtility<2>::PreReseed)
        .def("PostReseed", &MoveShallowWaterParticleUtility<2>::PostReseed)
        .def("ResetBoundaryConditions", &MoveShallowWaterParticleUtility<2>::ResetBoundaryConditions)
        .def("TransferLagrangianToEulerian",&MoveShallowWaterParticleUtility<2>::TransferLagrangianToEulerian)
        .def("CalculateVelOverElemSize", &MoveShallowWaterParticleUtility<2>::CalculateVelOverElemSize)
        .def("CalculateDeltaVariables", &MoveShallowWaterParticleUtility<2>::CalculateDeltaVariables)
        .def("CopyScalarVarToPreviousTimeStep", &MoveShallowWaterParticleUtility<2>::CopyScalarVarToPreviousTimeStep)
        .def("CopyVectorVarToPreviousTimeStep", &MoveShallowWaterParticleUtility<2>::CopyVectorVarToPreviousTimeStep)
        .def("ExecuteParticlesPrintingTool", &MoveShallowWaterParticleUtility<2>::ExecuteParticlesPrintingTool)
        ;

    py::class_< ShallowWaterUtilities > (m, "ShallowWaterUtilities")
        .def(py::init<>())
        .def("ComputeFreeSurfaceElevation", &ShallowWaterUtilities::ComputeFreeSurfaceElevation)
        .def("ComputeHeightFromFreeSurface", &ShallowWaterUtilities::ComputeHeightFromFreeSurface)
        .def("ComputeVelocity", &ShallowWaterUtilities::ComputeVelocity)
        .def("ComputeMomentum", &ShallowWaterUtilities::ComputeMomentum)
        .def("ComputeEnergy", &ShallowWaterUtilities::ComputeEnergy)
        .def("ComputeAccelerations", &ShallowWaterUtilities::ComputeAccelerations)
        .def("FlipScalarVariable", &ShallowWaterUtilities::FlipScalarVariable)
        .def("IdentifySolidBoundary", &ShallowWaterUtilities::IdentifySolidBoundary)
        .def("IdentifyWetDomain", &ShallowWaterUtilities::IdentifyWetDomain)
        .def("ResetDryDomain", &ShallowWaterUtilities::ResetDryDomain)
        .def("CopyFlag", &ShallowWaterUtilities::CopyFlag<ModelPart::NodesContainerType>)
        .def("CopyFlag", &ShallowWaterUtilities::CopyFlag<ModelPart::ElementsContainerType>)
        .def("CopyFlag", &ShallowWaterUtilities::CopyFlag<ModelPart::ConditionsContainerType>)
        .def("NormalizeVector", &ShallowWaterUtilities::NormalizeVector)
        .def("CopyVariableToPreviousTimeStep", &ShallowWaterUtilities::CopyVariableToPreviousTimeStep<Variable<double>&>)
        .def("CopyVariableToPreviousTimeStep", &ShallowWaterUtilities::CopyVariableToPreviousTimeStep<Variable<array_1d<double,3>>&>)
        .def("SetMinimumValue", &ShallowWaterUtilities::SetMinimumValue)
        .def("SetMeshZCoordinateToZero", &ShallowWaterUtilities::SetMeshZCoordinateToZero)
        .def("SetMeshZCoordinate", &ShallowWaterUtilities::SetMeshZCoordinate)
        ;

    py::class_< EstimateTimeStepUtility > (m, "EstimateTimeStepUtility")
        .def(py::init<ModelPart&, Parameters>())
        .def("Execute", &EstimateTimeStepUtility::Execute)
        ;

    py::class_< ReplicateModelPartUtility > (m, "ReplicateModelPartUtility")
        .def(py::init<ModelPart&, ModelPart&>())
        .def(py::init<ModelPart&, ModelPart&, bool>())
        .def("Replicate", &ReplicateModelPartUtility::Replicate)
        .def("TransferVariable", &ReplicateModelPartUtility::TransferVariable<Variable<double>>)
        .def("TransferVariable", &ReplicateModelPartUtility::TransferVariable<Variable<array_1d<double, 3>>>)
        .def("TransferNonHistoricalVariable", &ReplicateModelPartUtility::TransferNonHistoricalVariable<Variable<double>>)
        .def("TransferNonHistoricalVariable", &ReplicateModelPartUtility::TransferNonHistoricalVariable<Variable<array_1d<double, 3>>>)
        ;

    py::class_< PostProcessUtilities > (m, "PostProcessUtilities")
        .def(py::init<ModelPart&>())
        .def("DefineAuxiliaryProperties", &PostProcessUtilities::DefineAuxiliaryProperties)
        .def("AssignDryWetProperties", &PostProcessUtilities::AssignDryWetProperties)
        .def("RestoreDryWetProperties", &PostProcessUtilities::RestoreDryWetProperties)
        ;

    py::class_< BFECCConvectionUtility<2> > (m, "BFECCConvectionUtility")
        .def(py::init<ModelPart&>())
        .def(py::init<ModelPart&, Parameters>())
        .def("Convect", &BFECCConvectionUtility<2>::Convect<Variable<double>,double>)
        .def("Convect", &BFECCConvectionUtility<2>::Convect<Variable<array_1d<double,3>>,array_1d<double,3>>)
        .def("UpdateSearchDatabase", &BFECCConvectionUtility<2>::UpdateSearchDatabase)
        .def("ResetBoundaryConditions", &BFECCConvectionUtility<2>::ResetBoundaryConditions<Variable<double>>)
        .def("CopyVariableToPreviousTimeStep", &BFECCConvectionUtility<2>::CopyVariableToPreviousTimeStep<Variable<double>>)
        .def("CopyVariableToPreviousTimeStep", &BFECCConvectionUtility<2>::CopyVariableToPreviousTimeStep<Variable<array_1d<double,3>>>)
        ;

  }

}  // namespace Python.

} // Namespace Kratos
