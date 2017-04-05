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
#include "custom_utilities/move_particle_utility_pfem2.h" 

#include "custom_utilities/visualization.h" 
#include "custom_utilities/calculate_water_fraction.h"
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "custom_utilities/add_fixed_velocity_condition.h"

#include "custom_utilities/particle_utilities.h"
#include "custom_utilities/ulf_apply_bc_process.h"
#include "custom_utilities/pfem2_utilities.h"
#include "custom_utilities/mark_outer_nodes_process.h"
#include "custom_utilities/mark_fluid_process.h"
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


	class_< MoveParticleUtilityPFEM2<2> > ("MoveParticleUtilityPFEM22D", init<ModelPart& , int >())
	  .def("MountBin", &MoveParticleUtilityPFEM2<2>::MountBin)
	  .def("MoveParticles", &MoveParticleUtilityPFEM2<2>::MoveParticles)
	  .def("AccelerateParticlesWithoutMovingUsingDeltaVelocity", &MoveParticleUtilityPFEM2<2>::AccelerateParticlesWithoutMovingUsingDeltaVelocity)
	  .def("PreReseed", &MoveParticleUtilityPFEM2<2>::PreReseed)
	  .def("PostReseed", &MoveParticleUtilityPFEM2<2>::PostReseed)
	  .def("ResetBoundaryConditions", &MoveParticleUtilityPFEM2<2>::ResetBoundaryConditions)
	  .def("ResetBoundaryConditionsSlip", &MoveParticleUtilityPFEM2<2>::ResetBoundaryConditionsSlip)
	  .def("TransferLagrangianToEulerian",&MoveParticleUtilityPFEM2<2>::TransferLagrangianToEulerian)
	  .def("CalculateVelOverElemSize", &MoveParticleUtilityPFEM2<2>::CalculateVelOverElemSize)
	  .def("CalculateDeltaVelocity", &MoveParticleUtilityPFEM2<2>::CalculateDeltaVelocity)
	  .def("CopyVectorVarToPreviousTimeStep", &MoveParticleUtilityPFEM2<2>::CopyVectorVarToPreviousTimeStep)
	  .def("IntializeTransferTool", &MoveParticleUtilityPFEM2<2>::IntializeTransferTool)
	  .def("PreReseedUsingTopographicDomain", &MoveParticleUtilityPFEM2<2>::PreReseedUsingTopographicDomain)
	  .def("ExecuteParticlesPritingTool", &MoveParticleUtilityPFEM2<2>::ExecuteParticlesPritingTool)
	  .def("ExecuteParticlesPritingToolForDroppletsOnly", &MoveParticleUtilityPFEM2<2>::ExecuteParticlesPritingToolForDroppletsOnly)
	  .def("RotateParticlesAndDomainVelocities", &MoveParticleUtilityPFEM2<2>::RotateParticlesAndDomainVelocities)
	  ;    

	class_< MoveParticleUtilityPFEM2<3> > ("MoveParticleUtilityPFEM23D", init<ModelPart& , int >())
	  .def("MountBin", &MoveParticleUtilityPFEM2<3>::MountBin)
	  .def("MoveParticles", &MoveParticleUtilityPFEM2<3>::MoveParticles)
	  .def("AccelerateParticlesWithoutMovingUsingDeltaVelocity", &MoveParticleUtilityPFEM2<3>::AccelerateParticlesWithoutMovingUsingDeltaVelocity)
	  .def("PreReseed", &MoveParticleUtilityPFEM2<3>::PreReseed)
	  .def("PostReseed", &MoveParticleUtilityPFEM2<3>::PostReseed)
	  .def("ResetBoundaryConditions", &MoveParticleUtilityPFEM2<3>::ResetBoundaryConditions)
	  .def("ResetBoundaryConditionsSlip", &MoveParticleUtilityPFEM2<3>::ResetBoundaryConditionsSlip)
	  .def("TransferLagrangianToEulerian",&MoveParticleUtilityPFEM2<3>::TransferLagrangianToEulerian)
	  .def("CalculateVelOverElemSize", &MoveParticleUtilityPFEM2<3>::CalculateVelOverElemSize)
	  .def("CalculateDeltaVelocity", &MoveParticleUtilityPFEM2<3>::CalculateDeltaVelocity)
	  .def("CopyVectorVarToPreviousTimeStep", &MoveParticleUtilityPFEM2<3>::CopyVectorVarToPreviousTimeStep)
	  .def("IntializeTransferTool", &MoveParticleUtilityPFEM2<3>::IntializeTransferTool)
	  .def("PreReseedUsingTopographicDomain", &MoveParticleUtilityPFEM2<3>::PreReseedUsingTopographicDomain)
	  .def("ExecuteParticlesPritingTool", &MoveParticleUtilityPFEM2<3>::ExecuteParticlesPritingTool)
	  .def("ExecuteParticlesPritingToolForDroppletsOnly", &MoveParticleUtilityPFEM2<3>::ExecuteParticlesPritingToolForDroppletsOnly)
	  .def("AssignNodalVelocityUsingInletConditions", &MoveParticleUtilityPFEM2<3>::AssignNodalVelocityUsingInletConditions)
	  .def("RotateParticlesAndDomainVelocities", &MoveParticleUtilityPFEM2<3>::RotateParticlesAndDomainVelocities)
	  ;    
                    
      	class_<AddFixedVelocityCondition2D > ("AddFixedVelocityCondition2D", init<ModelPart& >())
	  .def("AddThem", &AddFixedVelocityCondition2D::AddThem)
	  ;
                    
        class_<AddWaterFixedVelocityCondition2D > ("AddWaterFixedVelocityCondition2D", init<ModelPart& >())
	  .def("AddThem", &AddWaterFixedVelocityCondition2D::AddThem)
	  ;
                    
        class_<AddFixedVelocityCondition3D > ("AddFixedVelocityCondition3D", init<ModelPart& >())
	  .def("AddThem", &AddFixedVelocityCondition3D::AddThem)
	  ;
                    
        class_<AddFixedPressureCondition2D > ("AddFixedPressureCondition2D", init<ModelPart& >())
	  .def("AddThem", &AddFixedPressureCondition2D::AddThem)
	  ;
                    
        class_<AddFixedPressureCondition3D > ("AddFixedPressureCondition3D", init<ModelPart& >())
	  .def("AddThem", &AddFixedPressureCondition3D::AddThem)
	  ;            
                    
	class_<VisualizationUtilities > ("VisualizationUtilities", init<>())
	  .def("VisualizationModelPart",&VisualizationUtilities::VisualizationModelPart)
	  ;
                    
        class_<CalculateWaterFraction<2> > ("CalculateWaterFraction2D", init<ModelPart& >())
	  .def("Calculate",&CalculateWaterFraction<2>::Calculate)
	  .def("CalculateWaterHeight",&CalculateWaterFraction<2>::CalculateWaterHeight)
	  .def("CalculateMeanCourant",&CalculateWaterFraction<2>::CalculateMeanCourant)
	  .def("CalculateMaxCourant",&CalculateWaterFraction<2>::CalculateMaxCourant)
	  .def("CalculateMaxCourantInNegativeElements",&CalculateWaterFraction<2>::CalculateMaxCourantInNegativeElements)
	  .def("CalculateForce",&CalculateWaterFraction<2>::CalculateForce)
	  ;
       
        class_<CalculateWaterFraction<3> > ("CalculateWaterFraction3D", init<ModelPart& >())
	  .def("Calculate",&CalculateWaterFraction<3>::Calculate)
	  .def("CalculateWaterHeight",&CalculateWaterFraction<3>::CalculateWaterHeight)
	  .def("CalculateMeanCourant",&CalculateWaterFraction<3>::CalculateMeanCourant)
	  .def("CalculateMaxCourant",&CalculateWaterFraction<3>::CalculateMaxCourant)
	  .def("CalculateMaxCourantInNegativeElements",&CalculateWaterFraction<3>::CalculateMaxCourantInNegativeElements)
	  .def("CalculateForce",&CalculateWaterFraction<3>::CalculateForce)
	  ;
	
	class_<ParticleUtils < 2 > >("ParticleUtils2D", init<>())
	  .def("VisualizationModelPart", &ParticleUtils < 2 > ::VisualizationModelPart)
	  .def("TransferToEulerianMesh", &ParticleUtils < 2 > ::TransferToEulerianMesh)
	  .def("MoveMesh_Streamlines_freesurfaceflows", &ParticleUtils < 2 > ::MoveMesh_Streamlines_freesurfaceflows)
	  .def("MoveLonelyNodes", &ParticleUtils < 2 > ::MoveLonelyNodes)	
	  .def("TransferToEulerianMeshShapeBased_aux", &ParticleUtils < 2 > ::TransferToEulerianMeshShapeBased_aux)
	  .def("TransferToEulerianMesh_Face_Heat_Flux", &ParticleUtils < 2 > ::TransferToEulerianMesh_Face_Heat_Flux)
	  ;
	
	
	class_<ParticleUtils < 3 > >("ParticleUtils3D", init<>())
	  .def("MoveMesh_Streamlines_freesurfaceflows", &ParticleUtils < 3 > ::MoveMesh_Streamlines_freesurfaceflows)
	  .def("TransferToEulerianMeshShapeBased_aux_3D", &ParticleUtils < 3 > ::TransferToEulerianMeshShapeBased_aux_3D)
	  .def("TransferToEulerianMesh_Face_Heat_Flux", &ParticleUtils < 3 > ::TransferToEulerianMesh_Face_Heat_Flux)	
	  .def("TransferToEulerianMesh", &ParticleUtils < 3 > ::TransferToEulerianMesh)	
	  .def("TransferToEulerianMesh_Face_Heat_Flux", &ParticleUtils < 3 > ::TransferToEulerianMesh_Face_Heat_Flux)	
	  .def("CalculateNormal", &ParticleUtils < 3 > ::CalculateNormal)	
	  ;
	
	class_<Pfem2ApplyBCProcess, bases<Process> >("Pfem2ApplyBCProcess", init<ModelPart&>());
	
	class_<Pfem2Utils>("Pfem2Utils", init<>())
	  .def("ApplyBoundaryConditions",&Pfem2Utils::ApplyBoundaryConditions)
	  .def("MarkOuterNodes",&Pfem2Utils::MarkOuterNodes)
	  .def("MoveLonelyNodes",&Pfem2Utils::MoveLonelyNodes)
	  .def("MarkExcessivelyCloseNodes",&Pfem2Utils::MarkExcessivelyCloseNodes)
	  .def("MarkNodesCloseToWall", &Pfem2Utils::MarkNodesCloseToWall)
	  .def("MarkNodesTouchingWall", &Pfem2Utils::MarkNodesTouchingWall)
	  .def("MarkNodesCloseToWallForBladder", &Pfem2Utils::MarkNodesCloseToWallForBladder)
	  .def("MarkNodesCloseToFS", &Pfem2Utils::MarkNodesCloseToFS)
	  .def ("MarkLonelyNodesForErasing", &Pfem2Utils::MarkLonelyNodesForErasing)
	  .def ("SaveReducedPart", &Pfem2Utils::SaveReducedPart)
	  ;
	class_<MarkOuterNodesProcess, bases<Process> >("MarkOuterNodesProcess",init<ModelPart&>())
	  .def("MarkOuterNodes",&MarkOuterNodesProcess::MarkOuterNodes)

	  ;
	
 	class_<MarkFluidProcess, bases<Process> >("MarkFluidProcess", init<ModelPart&>());
	
  }

}  // namespace Python.

} // Namespace Kratos

