/*
==============================================================================
KratosIncompressibleFluidApplication 
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
//   Date:                $Date: 2008-10-13 08:17:41 $
//   Revision:            $Revision: 1.20 $
//
//


// System includes 

// External includes 
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/calculate_forces_utilities.h" 
#include "custom_utilities/level_set_utilities.h" 
#include "custom_utilities/level_set_utilities_implicitextrapolation.h" 
#include "custom_utilities/edge_data.h"
#include "custom_utilities/fluid_solver.h" 
#include "custom_utilities/levelset_fluid_solver.h" 
#include "custom_utilities/Turbolence_Smagorinsky.h" 
#include "custom_utilities/pure_convection_edgebased.h" 
#include "custom_utilities/coupled_eulerian_ulf_utilities.h" 
#include "custom_utilities/elembased_distance_utilities.h" 
#include "custom_utilities/elembased_extrapolation_utilities.h" 
#include "custom_utilities/elembased_BC_utilities.h" 
#include "custom_utilities/assign_no_slip_condition.h" 
#include "custom_utilities/mark_for_refinement.h" 
#include "custom_utilities/parallel_extrapolation_utilities.h" 
#include "custom_utilities/wave_generator.h"
#include "custom_utilities/estimate_dt_utilities.h"
#include "custom_utilities/lagrangian_particle_utilities.h"


//#include "custom_utilities/edgebased_levelset.h"



#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"


namespace Kratos
{
	
namespace Python
{
  void  AddCustomUtilitiesToPython()
  {
	using namespace boost::python;

 class_<CoupledEulerianUlfUtils>("CoupledEulerianUlfUtils", init<>())
		.def("SavePseudoLagPart",&CoupledEulerianUlfUtils::SavePseudoLagPart)
		.def("ApplyProjDirichlet",&CoupledEulerianUlfUtils::ApplyProjDirichlet)
		.def("FindInterface",&CoupledEulerianUlfUtils::FindInterface)
		.def("FindIntersectionOfEdges",&CoupledEulerianUlfUtils::FindIntersectionOfEdges)
		.def("DisableSubdomain",&CoupledEulerianUlfUtils::DisableSubdomain)
		;

	  class_<CalculateForcesUtils>("CalculateForcesUtils", init<>())
		.def("CalculateForces3D",&CalculateForcesUtils::CalculateForces3D)
		.def("CalculateForces2D",&CalculateForcesUtils::CalculateForces2D)
		.def("CalculatePressureForces2D",&CalculateForcesUtils::CalculatePressureForces2D)
		;

	  class_<RefinementUtilities>("RefinementUtilities", init<>())
		.def("MarkForRefinement",&RefinementUtilities::MarkForRefinement)
                .def("SubscaleErrorEstimate",&RefinementUtilities::SubscaleErrorEstimate)
		;
		
	  class_<AssignNoSlipCondition>("AssignNoSlipCondition", init<>())
		.def("AssignNoSlipCondition2D",&AssignNoSlipCondition::AssignNoSlipCondition2D)
		;



	  class_<SmagorinskyTurbulentModel>("SmagorinskyTurbulentModel", init<>())
	    .def("CalculateTurbulentViscosity2D",&SmagorinskyTurbulentModel::CalculateTurbulentViscosity<2>)
	    .def("CalculateTurbulentViscosity3D",&SmagorinskyTurbulentModel::CalculateTurbulentViscosity<3>)
	    ;
	  
		
	  
	  class_<LevelSetUtilities>("LevelSetUtilities", init<>())
			  .def("RegenerateFluidModelPart",&LevelSetUtilities::RegenerateFluidModelPart)
			  .def("MarkNodesAsVisited",&LevelSetUtilities::MarkNodesAsVisited)
			  .def("SetDistanceToNegative",&LevelSetUtilities::SetDistanceToNegative)
			  .def("ExtrapolateVelocities",&LevelSetUtilities::ExtrapolateVelocities)  
			  .def("ExtrapolateVelocitiesByLayer",&LevelSetUtilities::ExtrapolateVelocitiesByLayer)  
			  .def("GenerateModelPart",&LevelSetUtilities::GenerateModelPart)  			  
			  .def("ImplicitExtrapolation_PreProcess",&LevelSetUtilities::ImplicitExtrapolation_PreProcess)  
			  .def("ImplicitExtrapolation_PostProcess",&LevelSetUtilities::ImplicitExtrapolation_PostProcess)  
			  .def("PrepareForInternalFluidDistanceComputation",&LevelSetUtilities::PrepareForInternalFluidDistanceComputation)  
			  .def("FluidDistanceComputation_FromBoundary",&LevelSetUtilities::FluidDistanceComputation_FromBoundary)
			  .def("ApplyMinimumExtrapolationPressureFix",&LevelSetUtilities::ApplyMinimumExtrapolationPressureFix)			  
			  ;
	  
	  class_<LevelSetUtilitiesImplicitExtrapolation>("LevelSetUtilitiesImplicitExtrapolation", init<>())
			  .def("RegenerateFluidModelPart",&LevelSetUtilitiesImplicitExtrapolation::RegenerateFluidModelPart)
			  .def("MarkNodesAsVisited",&LevelSetUtilitiesImplicitExtrapolation::MarkNodesAsVisited)
			  .def("SetDistanceToNegative",&LevelSetUtilitiesImplicitExtrapolation::SetDistanceToNegative)
			  .def("ExtrapolateVelocities",&LevelSetUtilitiesImplicitExtrapolation::ExtrapolateVelocities)  
			  .def("PrepareForInternalFluidDistanceComputation",&LevelSetUtilitiesImplicitExtrapolation::PrepareForInternalFluidDistanceComputation)  
			  .def("FluidDistanceComputation_FromBoundary",&LevelSetUtilitiesImplicitExtrapolation::FluidDistanceComputation_FromBoundary)
			  .def("ApplyFluidProperties",&LevelSetUtilitiesImplicitExtrapolation::ApplyFluidProperties)
			  ;

	  class_< ElemBasedDistanceUtilities>("ElemBasedDistanceUtilities", init<ModelPart& >() )
			  .def("IdentifyFreeSurface",&ElemBasedDistanceUtilities::IdentifyFreeSurface)
			  .def("MarkExternalAndMixedNodes",&ElemBasedDistanceUtilities::MarkExternalAndMixedNodes)
			  .def("MarkInternalAndMixedNodes",&ElemBasedDistanceUtilities::MarkInternalAndMixedNodes) 
			  .def("SaveScalarVariableToOldStep",&ElemBasedDistanceUtilities::SaveScalarVariableToOldStep)
			  .def("ChangeSignToDistance",&ElemBasedDistanceUtilities::ChangeSignToDistance)
			  .def("MarkNodesByDistance",&ElemBasedDistanceUtilities::MarkNodesByDistance)
			  ;

	  class_< ElemBasedExtrapolationUtilities > ("ElemBasedExtrapolationUtilities", init<ModelPart& >() )
			  .def("ExtrapolateVelocities",&ElemBasedExtrapolationUtilities::ExtrapolateVelocities)
			  ;
      
	  class_< ElemBasedBCUtilities >("ElemBasedBCUtilities", init<ModelPart& >() )
			.def("SetDividedElem_2D",&ElemBasedBCUtilities::SetDividedElem_2D)
			.def("SetPressureAndVelocityFixities",&ElemBasedBCUtilities::SetPressureAndVelocityFixities) 
			.def("FreePressureAndVelocity",&ElemBasedBCUtilities::FreePressureAndVelocity)
			.def("SetToZeroPressureAndVelocity",&ElemBasedBCUtilities::SetToZeroPressureAndVelocity)
			;

          typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
          typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
          typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
          // 	class_< ElementBasedNavierStokesSolver< 2, SparseSpaceType, LinearSolverType>,  boost::noncopyable >       ("ElementBasedNavierStokesSolver2D", init<ModelPart&>() )
//                           .def("ConstructSystemStructure",&ElementBasedNavierStokesSolver< 2, SparseSpaceType, LinearSolverType>::ConstructSystemStructure)
//                           .def("Clear",&ElementBasedNavierStokesSolver< 2, SparseSpaceType, LinearSolverType>::Clear)
//                           .def("SolveStep1",&ElementBasedNavierStokesSolver< 2, SparseSpaceType, LinearSolverType>::SolveStep1)
//                           .def("SolveStep2",&ElementBasedNavierStokesSolver< 2, SparseSpaceType, LinearSolverType>::SolveStep2)
//                           .def("SolveStep3",&ElementBasedNavierStokesSolver< 2, SparseSpaceType, LinearSolverType>::SolveStep3)
//                         ;
          class_< MatrixContainer< 2, SparseSpaceType>,  boost::noncopyable >   ("MatrixContainer2D", init<     >() )
                          .def("ConstructCSRVector",&MatrixContainer< 2, SparseSpaceType >::ConstructCSRVector)
                          .def("BuildCSRData",&MatrixContainer< 2, SparseSpaceType >::BuildCSRData)
                          .def("Clear",&MatrixContainer< 2, SparseSpaceType >::Clear)
                        ;

          class_< MatrixContainer< 3, SparseSpaceType>,  boost::noncopyable >   ("MatrixContainer3D", init<     >() )
                          .def("ConstructCSRVector",&MatrixContainer< 3, SparseSpaceType >::ConstructCSRVector)
                          .def("BuildCSRData",&MatrixContainer< 3, SparseSpaceType >::BuildCSRData)
                          .def("Clear",&MatrixContainer< 3, SparseSpaceType >::Clear)
                        ;

			
	  class_< FluidSolver< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType>,  boost::noncopyable >       ("FluidSolver2D", init<MatrixContainer< 2, SparseSpaceType>&, ModelPart& ,const double ,const double,const Vector, bool,double ,double ,double ,double ,bool  >() )
                          .def("Initialize",&FluidSolver< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::Initialize)
                          .def("ComputeTimeStep",&FluidSolver< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ComputeTimeStep)
                          .def("SolveStep1",&FluidSolver< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SolveStep1)
                          .def("SolveStep2",&FluidSolver< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SolveStep2)
                          .def("SolveStep3",&FluidSolver< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SolveStep3)
                          .def("ComputeTimeStep",&FluidSolver< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ComputeTimeStep)
                          .def("CalculateNormals",&FluidSolver< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::CalculateNormals)
                          .def("UpdateFixedVelocityValues",&FluidSolver< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::UpdateFixedVelocityValues)
                        ;

	  class_< FluidSolver< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType>,  boost::noncopyable >       ("FluidSolver3D", init<MatrixContainer< 3, SparseSpaceType>&, ModelPart& ,const double ,const double,const Vector, bool,double ,double ,double ,double ,bool  >() )
                          .def("Initialize",&FluidSolver< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::Initialize)
                          .def("ComputeTimeStep",&FluidSolver< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ComputeTimeStep)
                          .def("SolveStep1",&FluidSolver< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SolveStep1)
                          .def("SolveStep2",&FluidSolver< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SolveStep2)
                          .def("SolveStep3",&FluidSolver< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SolveStep3)
                          .def("ComputeTimeStep",&FluidSolver< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ComputeTimeStep)
                          .def("CalculateNormals",&FluidSolver< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::CalculateNormals)
                          .def("UpdateFixedVelocityValues",&FluidSolver< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::UpdateFixedVelocityValues)
			  ;
	  
	  
	  class_< LevelSetFluidSolver< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType>,  boost::noncopyable >       ("LevelSetFluidSolver2D", init<MatrixContainer< 2, SparseSpaceType>&, ModelPart&, bool, bool >() )
			  .def("Initialize",&LevelSetFluidSolver< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::Initialize)
			  .def("ComputeTimeStep",&LevelSetFluidSolver< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ComputeTimeStep)
			  .def("SolveStep1",&LevelSetFluidSolver< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SolveStep1)
			  .def("SolveStep2",&LevelSetFluidSolver< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SolveStep2)
			  .def("SolveStep3",&LevelSetFluidSolver< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SolveStep3)
			  .def("SolveStep4",&LevelSetFluidSolver< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SolveStep4)
			  .def("ExtrapolateVelocities",&LevelSetFluidSolver< 2,MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ExtrapolateVelocities)
			  .def("MarkExternalAndMixedNodes",&LevelSetFluidSolver< 2,MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::MarkExternalAndMixedNodes)
			  .def("MarkInternalAndMixedNodes",&LevelSetFluidSolver< 2,MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::MarkInternalAndMixedNodes) 
			  .def("SaveScalarVariableToOldStep",&LevelSetFluidSolver< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SaveScalarVariableToOldStep)
			  .def("ChangeSignToDistance",&LevelSetFluidSolver< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ChangeSignToDistance)
			  .def("MarkNodesByDistance",&LevelSetFluidSolver< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::MarkNodesByDistance)
			  .def("CalculateForces",&LevelSetFluidSolver< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::CalculateForces)
			  .def("CalculateVariablesDistribution",&LevelSetFluidSolver< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::CalculateVariablesDistribution)
        		  ;

	  class_< LevelSetFluidSolver< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType>,  boost::noncopyable >       ("LevelSetFluidSolver3D", init<MatrixContainer< 3, SparseSpaceType>&, ModelPart&, bool, bool >() )
			  .def("Initialize",&LevelSetFluidSolver< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::Initialize)
			  .def("ComputeTimeStep",&LevelSetFluidSolver< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ComputeTimeStep)
			  .def("SolveStep1",&LevelSetFluidSolver< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SolveStep1)
			  .def("SolveStep2",&LevelSetFluidSolver< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SolveStep2)
			  .def("SolveStep3",&LevelSetFluidSolver< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SolveStep3)
			  .def("SolveStep4",&LevelSetFluidSolver< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SolveStep4)
			  .def("ExtrapolateVelocities",&LevelSetFluidSolver< 3,MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ExtrapolateVelocities)
			  .def("MarkExternalAndMixedNodes",&LevelSetFluidSolver< 3,MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::MarkExternalAndMixedNodes)
			  .def("MarkInternalAndMixedNodes",&LevelSetFluidSolver< 3,MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::MarkInternalAndMixedNodes) 
			  .def("SaveScalarVariableToOldStep",&LevelSetFluidSolver< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SaveScalarVariableToOldStep)
			  .def("ChangeSignToDistance",&LevelSetFluidSolver< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ChangeSignToDistance)
			  .def("MarkNodesByDistance",&LevelSetFluidSolver< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::MarkNodesByDistance)
			  .def("CalculateForces",&LevelSetFluidSolver< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::CalculateForces)
			  .def("CalculateVariablesDistribution",&LevelSetFluidSolver< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::CalculateVariablesDistribution)
        		;


	 class_< PureConvectionEdgeBased< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType>,  boost::noncopyable >       ("PureConvectionEdgeBased2D", init<MatrixContainer< 2, SparseSpaceType>&, ModelPart& >() )
                          .def("Initialize",&PureConvectionEdgeBased< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::Initialize)
                          .def("Solve",&PureConvectionEdgeBased< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::Solve)
			  .def("ComputeTimeStep",&PureConvectionEdgeBased< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ComputeTimeStep)
                        ;


	 class_< PureConvectionEdgeBased< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType>,  boost::noncopyable >       ("PureConvectionEdgeBased3D", init<MatrixContainer< 3, SparseSpaceType>&, ModelPart& >() )
                          .def("Initialize",&PureConvectionEdgeBased< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::Initialize)
                          .def("Solve",&PureConvectionEdgeBased< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::Solve)
			  .def("ComputeTimeStep",&PureConvectionEdgeBased< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ComputeTimeStep)
                        ;

	class_< ParallelExtrapolationUtilities<  3 >,  boost::noncopyable >       ("ParallelExtrapolationUtilities3D", init<>() )
                          .def("ExtrapolateVelocity",&ParallelExtrapolationUtilities< 3 >::ExtrapolateVelocity)
                          .def("ExtrapolatePressureProjection",&ParallelExtrapolationUtilities< 3 >::ExtrapolatePressureProjection)
                          .def("AssignFreeSurfacePressure",&ParallelExtrapolationUtilities< 3 >::AssignFreeSurfacePressure)
                        ;

	class_< WaveGenerator,  boost::noncopyable >       ("WaveGenerator", init<>() )
                          .def("GenerateWaveXYPlane",&WaveGenerator::GenerateWaveXYPlane)
                        ;

        class_< EstimateDtUtil< 2 >, boost::noncopyable > ("EstimateDt2D", init<ModelPart&>() )
                          .def("EstimateDt", &EstimateDtUtil<2>::EstimateDt)
                        ;

        class_< EstimateDtUtil< 3 >, boost::noncopyable > ("EstimateDt3D", init<ModelPart&>() )
                          .def("EstimateDt", &EstimateDtUtil<3>::EstimateDt)
                        ;
			
// 	class_< ElementBasedNavierStokesSolver< 2, SparseSpaceType, LinearSolverType>,  boost::noncopyable >       ("ElementBasedNavierStokesSolver2D", init<ModelPart&>() )
//                           .def("ConstructSystemStructure",&ElementBasedNavierStokesSolver< 2, SparseSpaceType, LinearSolverType>::ConstructSystemStructure)
//                           .def("Clear",&ElementBasedNavierStokesSolver< 2, SparseSpaceType, LinearSolverType>::Clear)
//                           .def("SolveStep1",&ElementBasedNavierStokesSolver< 2, SparseSpaceType, LinearSolverType>::SolveStep1)
//                           .def("SolveStep2",&ElementBasedNavierStokesSolver< 2, SparseSpaceType, LinearSolverType>::SolveStep2)
//                           .def("SolveStep3",&ElementBasedNavierStokesSolver< 2, SparseSpaceType, LinearSolverType>::SolveStep3)
// 			  .def("CalculateProjection",&ElementBasedNavierStokesSolver< 2, SparseSpaceType, LinearSolverType>::CalculateProjection)
// 			;

//	  class_< EdgeBasedLevelSet< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType>,  boost::noncopyable >       ("EdgeBasedLevelSet2D", init< MatrixContainer< 2, SparseSpaceType>&, ModelPart& >() )
//			  .def("Initialize",&EdgeBasedLevelSet< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::Initialize)
//			  .def("ComputeTimeStep",&EdgeBasedLevelSet< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ComputeTimeStep)
//			  .def("SolveStep1",&EdgeBasedLevelSet< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SolveStep1)
//			  .def("SolveStep2",&EdgeBasedLevelSet< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SolveStep2)
//			  .def("SolveStep3",&EdgeBasedLevelSet< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SolveStep3)
//			  .def("ExtrapolateValues",&EdgeBasedLevelSet< 2,MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ExtrapolateValues)
//			  .def("MarkExternalAndMixedNodes",&EdgeBasedLevelSet< 2,MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::MarkExternalAndMixedNodes)
//			  .def("MarkInternalAndMixedNodes",&EdgeBasedLevelSet< 2,MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::MarkInternalAndMixedNodes)
//			  .def("SaveScalarVariableToOldStep",&EdgeBasedLevelSet< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SaveScalarVariableToOldStep)
//			  .def("ChangeSignToDistance",&EdgeBasedLevelSet< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ChangeSignToDistance)
//			  .def("MarkNodesByDistance",&EdgeBasedLevelSet< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::MarkNodesByDistance)
////			  .def("CalculateForces",&EdgeBasedLevelSet< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::CalculateForces)
//			  ;

        class_<LagrangianParticleUtils < 2 > >("LagrangianUtils2D", init<>())
                    .def("StreamlineMove", &LagrangianParticleUtils < 2 > ::StreamlineMove)
                    .def("Reseed", &LagrangianParticleUtils < 2 > ::Reseed)
                    .def("VisualizationModelPart", &LagrangianParticleUtils < 2 > ::VisualizationModelPart)
                    .def("TransferToEulerianMesh", &LagrangianParticleUtils < 2 > ::TransferToEulerianMesh)
                    .def("RestartStep", &LagrangianParticleUtils < 2 > ::RestartStep)
                    ;

  }
	
}  // namespace Python.

} // Namespace Kratos

