//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Kazem Kamran
//                   Jordi Rubio
//                   Riccardo Rossi
//

// System includes

// External includes
#include "pybind11/pybind11.h"


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/assign_environment_condition.h"
#include "custom_utilities/estimate_time_step.h"
#include "custom_utilities/particle_levelset_utilities.h"
#include "custom_utilities/biphasic_filling_utilities.h"
#include "includes/deprecated_variables.h"
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"



namespace Kratos
{

namespace Python
{


void  AddCustomUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;


//     typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
//     typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
//     typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;


    py::class_<AssignEnvironmentCondition > (m, "AssignEnvironmentCondition")
    .def(py::init<>())
    .def("AssignCondition", &AssignEnvironmentCondition::AssignCondition)
    ;
    py::class_<EstimateTimeStep < 3 >, typename EstimateTimeStep < 3 >::Pointer > (m, "EstimateTimeStep3D")
    .def(py::init<>())
    .def("ComputeDt", &EstimateTimeStep < 3 >::ComputeDt)
    .def("ComputeSolidificationCoolingDt", &EstimateTimeStep < 3 >::ComputeSolidificationCoolingDt)
    .def("EstimateSolidificationTime", &EstimateTimeStep < 3 >::EstimateSolidificationTime)
	.def("EstimateSolidificationTimeNoVirtualMould", &EstimateTimeStep < 3 >::EstimateSolidificationTimeNoVirtualMould)
    .def("CheckStopTemperature", &EstimateTimeStep < 3 >::CheckStopTemperature)
    .def("ComputeSurfaceWaveDt", &EstimateTimeStep < 3 >::ComputeSurfaceWaveDt)
    .def("CheckIsInTransition", &EstimateTimeStep < 3 >::CheckIsInTransition)
	.def("EstimateCoolingTime", &EstimateTimeStep < 3 >::EstimateCoolingTime)
	.def("CheckMinTemperature", &EstimateTimeStep < 3 >::CheckMinTemperature)
    ;

    py::class_<ParticleLevelSetUtils < 2 >, typename ParticleLevelSetUtils < 2 >::Pointer >(m, "ParticleLevelSetUtils2D")
    .def(py::init<>())
    .def("Seed", &ParticleLevelSetUtils < 2 > ::Seed)
    .def("StreamlineMove", &ParticleLevelSetUtils < 2 > ::StreamlineMove)
    .def("VisualizationModelPart", &ParticleLevelSetUtils < 2 > ::VisualizationModelPart)
    .def("FindMaxMinEdgeSize", &ParticleLevelSetUtils < 2 > ::FindMaxMinEdgeSize)
    .def("ResetParticleRadius", &ParticleLevelSetUtils < 2 > ::ResetParticleRadius)
    .def("ParticleLevelSetCorrection", &ParticleLevelSetUtils < 2 > ::ParticleLevelSetCorrection)
    .def("ParticleReseeding", &ParticleLevelSetUtils < 2 > ::ParticleReseeding)
    ;

    py::class_<ParticleLevelSetUtils < 3 >, typename ParticleLevelSetUtils < 3 >::Pointer >(m, "ParticleLevelSetUtils3D")
    .def(py::init<>())
    .def("Seed", &ParticleLevelSetUtils < 3 > ::Seed)
    .def("StreamlineMove", &ParticleLevelSetUtils < 3 > ::StreamlineMove)
    .def("VisualizationModelPart", &ParticleLevelSetUtils < 3 > ::VisualizationModelPart)
    .def("FindMaxMinEdgeSize", &ParticleLevelSetUtils < 3 > ::FindMaxMinEdgeSize)
    .def("ResetParticleRadius", &ParticleLevelSetUtils < 3 > ::ResetParticleRadius)
    .def("ParticleLevelSetCorrection", &ParticleLevelSetUtils < 3 > ::ParticleLevelSetCorrection)
    .def("ParticleReseeding", &ParticleLevelSetUtils < 3 > ::ParticleReseeding)
    ;
    py::class_<BiphasicFillingUtilities, BiphasicFillingUtilities::Pointer > (m, "BiphasicFillingUtilities")
    .def(py::init<>())
    .def("CreateAutoExitAssignAirSmagorinsky", &BiphasicFillingUtilities::CreateAutoExitAssignAirSmagorinsky)
	.def("AssignSmoothBoundaryAirExit", &BiphasicFillingUtilities::AssignSmoothBoundaryAirExit)
	.def("ApplyFluidProperties", &BiphasicFillingUtilities::ApplyFluidProperties)
	.def("DistanceFarRegionCorrection", &BiphasicFillingUtilities::DistanceFarRegionCorrection)
	.def("VolumeCorrection", &BiphasicFillingUtilities::VolumeCorrection)
    .def("ComputeFillPercentage", &BiphasicFillingUtilities::ComputeFillPercentage)
	.def("ComputeNetInletVolume", &BiphasicFillingUtilities::ComputeNetInletVolume)
	.def("ComputeNodalVolume", &BiphasicFillingUtilities::ComputeNodalVolume)
    .def("ApplyVelocityLimitation", &BiphasicFillingUtilities::ApplyVelocityLimitation)
	.def("LastStepExtrapolations", &BiphasicFillingUtilities::LastStepExtrapolations)
	.def("SolidificationDuringFilling", &BiphasicFillingUtilities::SolidificationDuringFilling)
	.def("ViscosityBasedSolidification", &BiphasicFillingUtilities::ViscosityBasedSolidification)
	//.def("MacroPorosityToShrinkageComputation", &BiphasicFillingUtilities::MacroPorosityToShrinkageComputation)
	.def("ComputePosetiveVolume", &BiphasicFillingUtilities::ComputePosetiveVolume)
	.def("PosetiveVolumeCorrection", &BiphasicFillingUtilities::PosetiveVolumeCorrection)
	.def("ApplyTemperatureLimitation", &BiphasicFillingUtilities::ApplyTemperatureLimitation)
	.def("CheckIfAllNodesAreWet", &BiphasicFillingUtilities::CheckIfAllNodesAreWet)
	.def("ComputeWetVolume", &BiphasicFillingUtilities::ComputeWetVolume)
	.def("ComputePartVolume", &BiphasicFillingUtilities::ComputePartVolume)
	.def("ComputePartArea", &BiphasicFillingUtilities::ComputePartArea)
	.def("ComputePartInletArea", &BiphasicFillingUtilities::ComputePartInletArea)
	.def("ComputePartMaxh", &BiphasicFillingUtilities::ComputePartMaxh)
	.def("ComputePartAvgh", &BiphasicFillingUtilities::ComputePartAvgh)
    ;


}





}  // namespace Python.

} // Namespace Kratos

