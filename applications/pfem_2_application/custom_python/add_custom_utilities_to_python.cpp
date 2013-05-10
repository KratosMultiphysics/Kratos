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
#include "custom_utilities/move_particle_utility_diff.h" 
#include "custom_utilities/visualization.h" 
#include "custom_utilities/calculate_water_fraction.h"
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "custom_utilities/add_fixed_velocity_condition.h"

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


		class_< MoveParticleUtilityDiff<2> > ("MoveParticleUtilityDiff2D", init<ModelPart& , int >())
                    .def("MountBinDiff", &MoveParticleUtilityDiff<2>::MountBinDiff)
                    .def("MoveParticlesDiff", &MoveParticleUtilityDiff<2>::MoveParticlesDiff)
                    .def("AccelerateParticlesWithoutMoving", &MoveParticleUtilityDiff<2>::AccelerateParticlesWithoutMoving)
                    .def("AccelerateParticlesWithoutMovingUsingDeltaVelocity", &MoveParticleUtilityDiff<2>::AccelerateParticlesWithoutMovingUsingDeltaVelocity)
                    .def("PreReseed", &MoveParticleUtilityDiff<2>::PreReseed)
                    .def("PostReseed", &MoveParticleUtilityDiff<2>::PostReseed)
                    .def("ResetBoundaryConditions", &MoveParticleUtilityDiff<2>::ResetBoundaryConditions)
                    .def("TransferLagrangianToEulerian",&MoveParticleUtilityDiff<2>::TransferLagrangianToEulerian)
                    .def("ReplaceParticlesVelocityAndDistance",&MoveParticleUtilityDiff<2>::ReplaceParticlesVelocityAndDistance)
                    .def("CalculateSurfaceNormal", &MoveParticleUtilityDiff<2>::CalculateSurfaceNormal)
                    .def("FlagSplittedElementsAndTheirNodes", &MoveParticleUtilityDiff<2>::FlagSplittedElementsAndTheirNodes)
                    .def("CalculateVelOverElemSize", &MoveParticleUtilityDiff<2>::CalculateVelOverElemSize)
                    .def("ApplyDistanceFunctionCorrection", &MoveParticleUtilityDiff<2>::ApplyDistanceFunctionCorrection)
                    .def("CalculateDeltaVelocity", &MoveParticleUtilityDiff<2>::CalculateDeltaVelocity)
                    ; 
                    
		class_< MoveParticleUtilityDiff<3> > ("MoveParticleUtilityDiff3D", init<ModelPart& , int >())
                    .def("MountBinDiff", &MoveParticleUtilityDiff<3>::MountBinDiff)
                    .def("MoveParticlesDiff", &MoveParticleUtilityDiff<3>::MoveParticlesDiff)
                    .def("AccelerateParticlesWithoutMoving", &MoveParticleUtilityDiff<3>::AccelerateParticlesWithoutMoving)
                    .def("AccelerateParticlesWithoutMovingUsingDeltaVelocity", &MoveParticleUtilityDiff<3>::AccelerateParticlesWithoutMovingUsingDeltaVelocity)
                    .def("PreReseed", &MoveParticleUtilityDiff<3>::PreReseed)
                    .def("PostReseed", &MoveParticleUtilityDiff<3>::PostReseed)
                    .def("ResetBoundaryConditions", &MoveParticleUtilityDiff<3>::ResetBoundaryConditions)
                    .def("TransferLagrangianToEulerian",&MoveParticleUtilityDiff<3>::TransferLagrangianToEulerian)
                    .def("ReplaceParticlesVelocityAndDistance",&MoveParticleUtilityDiff<3>::ReplaceParticlesVelocityAndDistance)
                    .def("CalculateSurfaceNormal", &MoveParticleUtilityDiff<3>::CalculateSurfaceNormal)
                    .def("FlagSplittedElementsAndTheirNodes", &MoveParticleUtilityDiff<3>::FlagSplittedElementsAndTheirNodes)
                    .def("CalculateVelOverElemSize", &MoveParticleUtilityDiff<3>::CalculateVelOverElemSize)
                    .def("ApplyDistanceFunctionCorrection", &MoveParticleUtilityDiff<3>::ApplyDistanceFunctionCorrection)
                    .def("CalculateDeltaVelocity", &MoveParticleUtilityDiff<3>::CalculateDeltaVelocity)
                    ;            
                    
                    
                    
      	class_<AddFixedVelocityCondition2D > ("AddFixedVelocityCondition2D", init<ModelPart& >())
                    .def("AddThem", &AddFixedVelocityCondition2D::AddThem)
                    ;
                    
        class_<AddFixedVelocityCondition3D > ("AddFixedVelocityCondition3D", init<ModelPart& >())
                    .def("AddThem", &AddFixedVelocityCondition3D::AddThem)
                    ;
                    
		class_<VisualizationUtilities2D > ("VisualizationUtilities2D", init<>())
                    .def("VisualizationModelPart",&VisualizationUtilities2D::VisualizationModelPart)
                    ;
                    
        class_<CalculateWaterFraction > ("CalculateWaterFraction", init<ModelPart& >())
					.def("Calculate",&CalculateWaterFraction::Calculate)
					.def("CalculateWaterHeight",&CalculateWaterFraction::CalculateWaterHeight)
					.def("CalculateMeanCourant",&CalculateWaterFraction::CalculateMeanCourant)
					.def("CalculateMaxCourant",&CalculateWaterFraction::CalculateMaxCourant)
					.def("CalculateForce",&CalculateWaterFraction::CalculateForce)
                    ;
                    

  }
	




}  // namespace Python.

} // Namespace Kratos

