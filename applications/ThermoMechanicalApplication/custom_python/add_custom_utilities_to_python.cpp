/*
==============================================================================
KratosTestApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/

//
//   Project Name:        Kratos
//   Last modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: 1.2 $
//
//


// System includes

// External includes
#include <boost/python.hpp>


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


void  AddCustomUtilitiesToPython()
{
    using namespace boost::python;


    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;

    
    class_<AssignEnvironmentCondition > ("AssignEnvironmentCondition", init<>())
    .def("AssignCondition", &AssignEnvironmentCondition::AssignCondition)
    ;    
    class_<EstimateTimeStep < 3 > > ("EstimateTimeStep3D", init<>())
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
    
    class_<ParticleLevelSetUtils < 2 > >("ParticleLevelSetUtils2D", init<>())
    .def("Seed", &ParticleLevelSetUtils < 2 > ::Seed)
    .def("StreamlineMove", &ParticleLevelSetUtils < 2 > ::StreamlineMove)
    .def("VisualizationModelPart", &ParticleLevelSetUtils < 2 > ::VisualizationModelPart)
    .def("FindMaxMinEdgeSize", &ParticleLevelSetUtils < 2 > ::FindMaxMinEdgeSize)    
    .def("ResetParticleRadius", &ParticleLevelSetUtils < 2 > ::ResetParticleRadius)
    .def("ParticleLevelSetCorrection", &ParticleLevelSetUtils < 2 > ::ParticleLevelSetCorrection) 
    .def("ParticleReseeding", &ParticleLevelSetUtils < 2 > ::ParticleReseeding)
    ; 
    
    class_<ParticleLevelSetUtils < 3 > >("ParticleLevelSetUtils3D", init<>())
    .def("Seed", &ParticleLevelSetUtils < 3 > ::Seed)
    .def("StreamlineMove", &ParticleLevelSetUtils < 3 > ::StreamlineMove)
    .def("VisualizationModelPart", &ParticleLevelSetUtils < 3 > ::VisualizationModelPart)
    .def("FindMaxMinEdgeSize", &ParticleLevelSetUtils < 3 > ::FindMaxMinEdgeSize)      
    .def("ResetParticleRadius", &ParticleLevelSetUtils < 3 > ::ResetParticleRadius)
    .def("ParticleLevelSetCorrection", &ParticleLevelSetUtils < 3 > ::ParticleLevelSetCorrection) 
    .def("ParticleReseeding", &ParticleLevelSetUtils < 3 > ::ParticleReseeding)
    ;    
    class_<BiphasicFillingUtilities > ("BiphasicFillingUtilities", init<>())
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
	.def("MacroPorosityToShrinkageComputation", &BiphasicFillingUtilities::MacroPorosityToShrinkageComputation) 
	.def("ComputePosetiveVolume", &BiphasicFillingUtilities::ComputePosetiveVolume)
	.def("PosetiveVolumeCorrection", &BiphasicFillingUtilities::PosetiveVolumeCorrection)
	.def("ApplyTemperatureLimitation", &BiphasicFillingUtilities::ApplyTemperatureLimitation) 
	.def("CheckIfAllNodesAreWet", &BiphasicFillingUtilities::CheckIfAllNodesAreWet) 
	.def("ComputeWetVolume", &BiphasicFillingUtilities::ComputeWetVolume)
	
    ; 


}





}  // namespace Python.

} // Namespace Kratos

