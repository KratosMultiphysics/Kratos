/*
==============================================================================
KratosPFEMApplication 
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
//   Last modified by:    $Author: rrossi $
//   Date:                $Date: 2008-03-05 09:39:14 $
//   Revision:            $Revision: 1.6 $
//
//


// System includes 

// External includes 
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/pfem_utilities.h"
#include "custom_utilities/normal_to_wall_calculation_utils.h"
#include "custom_utilities/volume_correction_utils.h" 
#include "custom_utilities/lagrangian_utilities.h" 
#include "custom_utilities/nist_utilities.h" 
#include "custom_utilities/erosion_utilities.h" 
#include "custom_utilities/drag_utilities.h" 
//#include "custom_utilities/streamline_utils.h" 

#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "custom_utilities/lagrangian_euler_solver.h"
#include "custom_utilities/exact_dt_estimate_utilities.h"
namespace Kratos
{
	
namespace Python
{
	/*void GenerateModelPart(NistUtils& NistUtils,ModelPart& origin_model_part,ModelPart& destination_model_part,unsigned int domain_size )
	{
		if(domain_size == 2)
		{
			NistUtils.GenerateModelPart(origin_model_part, destination_model_part,
				KratosComponents<Element>::Get("ConvDiff2D"),
				KratosComponents<Condition>::Get("ThermalFace2D")	); 
		}
		else if(domain_size == 3)
		{
			NistUtils.GenerateModelPart(origin_model_part, destination_model_part,
				KratosComponents<Element>::Get("ConvDiff3D"),
				KratosComponents<Condition>::Get("ThermalFace3D")	); 
		}
	}*/
	
  void  AddCustomUtilitiesToPython()
  {
	using namespace boost::python;

	  class_<PfemUtils>("PfemUtils", init<>())
		.def("ApplyBoundaryConditions",&PfemUtils::ApplyBoundaryConditions)
		.def("ApplyMinimalPressureConditions",&PfemUtils::ApplyMinimalPressureConditions)
		.def("EstimateDeltaTime",&PfemUtils::EstimateDeltaTime)
		.def("MarkOuterNodes",&PfemUtils::MarkOuterNodes)
		.def("CalculateSpatialALEAcceleration2D",&PfemUtils::CalculateSpatialALEAcceleration<2>)
		.def("CalculateSpatialALEAcceleration3D",&PfemUtils::CalculateSpatialALEAcceleration<3>)
		.def("QuasiLagrangianMove",&PfemUtils::QuasiLagrangianMove)
		.def("Predict",&PfemUtils::Predict)
		.def("MoveLonelyNodes",&PfemUtils::MoveLonelyNodes)
		.def("CalculateVolume",&PfemUtils::CalculateVolume)
		.def("CalculateNodalMass",&PfemUtils::CalculateNodalMass)
		.def("IdentifyFluidNodes",&PfemUtils::IdentifyFluidNodes)
		.def("ReduceTimeStep",&PfemUtils::ReduceTimeStep)
		.def("MarkExcessivelyCloseNodes",&PfemUtils::MarkExcessivelyCloseNodes)
		.def("MarkNodesTouchingWall",&PfemUtils::MarkNodesTouchingWall)
		.def("MarkNodesTouchingInterface",&PfemUtils::MarkNodesTouchingInterface)
		.def("InterfaceDetecting",&PfemUtils::InterfaceDetecting)
		.def("ChangeWallWaterFlag",&PfemUtils::ChangeWallWaterFlag)
		.def("ChangeInterfaceWaterFlag",&PfemUtils::ChangeInterfaceWaterFlag)
		.def("ColourAirWaterElement",&PfemUtils::ColourAirWaterElement)
		.def("AssignNearBoundaryH",&PfemUtils::AssignNearBoundaryH)
		.def("MoveNodes",&PfemUtils::MoveNodes)
		.def("AssignMeshVelocity",&PfemUtils::AssignMeshVelocity)
                .def("CFLdeltaT",&PfemUtils::CFLdeltaT)
                .def("ExplicitDeltaT",&PfemUtils::ExplicitDeltaT)
		.def("CheckInvertedElements",&PfemUtils::CheckInvertedElements)
                
		 ;
	  
	  class_<NormalToWallCalculationUtils>("NormalToWallCalculationUtils", init<>())
		.def("CalculateNormalToWall",&NormalToWallCalculationUtils::CalculateNormalToWall)
		;

	  class_<VolumeCorrectionUtils>("VolumeCorrectionUtils", init<>())
		.def("CorrectVolume",&VolumeCorrectionUtils::CorrectVolume)
		;
	
	  class_<LagrangianUtils>("LagrangianUtils", init<>())
		.def("ExplicitLagrangianPrediction",&LagrangianUtils::ExplicitLagrangianPrediction)
		.def("ImplicitLagrangianPrediction",&LagrangianUtils::ImplicitLagrangianPrediction)
		.def("CalculateStep1DisplacementCorrection",&LagrangianUtils::CalculateStep1DisplacementCorrection)
		.def("CalculateFinalDisplacementCorrection",&LagrangianUtils::CalculateFinalDisplacementCorrection)
		;

	  class_<NistUtils>("NistUtils", init<>())
		//.def("GenerateModelPart",GenerateModelPart)
		.def("ApplyInitialTemperature",&NistUtils::ApplyInitialTemperature)
		.def("FindFluidLevel",&NistUtils::FindFluidLevel)
		;

	  class_<ErosionUtils<2>,  boost::noncopyable>  ("ErosionUtils2D", init< >())
		.def("CheckErosionableNodes",&ErosionUtils<2>::CheckErosionableNodes)
		.def("SetErosionableNodes", &ErosionUtils<2>::SetErosionableNodes)
		;

	  class_<ErosionUtils<3>,  boost::noncopyable>  ("ErosionUtils3D", init< >())
		.def("CheckErosionableNodes",&ErosionUtils<3>::CheckErosionableNodes)
		.def("SetErosionableNodes", &ErosionUtils<3>::SetErosionableNodes)
		;

		typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
		typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
		typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;

        class_< LagrangianEulerSolver< 2, SparseSpaceType, LinearSolverType >,  boost::noncopyable >	("LagrangianEulerSolver2D", init<	>() )
               .def("Predict",&LagrangianEulerSolver< 2, SparseSpaceType, LinearSolverType >::Predict)
               .def("SolveConvectionStep",&LagrangianEulerSolver< 2, SparseSpaceType, LinearSolverType >::SolveConvectionStep)
               .def("ConstructLaplacianSystem",&LagrangianEulerSolver< 2, SparseSpaceType, LinearSolverType >::ConstructLaplacianSystem)
               .def("BuildAndSolveLaplacianSystem",&LagrangianEulerSolver< 2, SparseSpaceType, LinearSolverType >::BuildAndSolveLaplacianSystem)
               .def("ClearLaplacianSystem",&LagrangianEulerSolver< 2, SparseSpaceType, LinearSolverType >::ClearLaplacianSystem)
        ;


	  class_<ExactDtEstimateUtilities>("ExactDtEstimateUtilities", init<>())
		.def("CubicExactDt",&ExactDtEstimateUtilities::CubicExactDt)
		;

	  class_<DragUtils>("DragUtils", init<>())
		.def("CalculateFluidDrag",&DragUtils::CalculateFluidDrag)
		.def("AddDrag",&DragUtils::AddDrag)
		;

 //         class_<StreamlineUtils<2>,  boost::noncopyable>  ("StreamlineUtils2D", init< >())
//		.def("LagrangianMoveBackAndForward",&StreamlineUtils<2>::LagrangianMoveBackAndForward)
//		 ;


  }
	




}  // namespace Python.

} // Namespace Kratos

