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
//   Last modified by:    $Author: kazem $
//   Date:                $Date: 2008-12-15 10:10:27 $
//   Revision:            $Revision: 1.8 $
//
//


// System includes 


// External includes 
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/timer.hpp> 


// Project includes
#include "includes/define.h"
#include "custom_python/add_custom_strategies_to_python.h"

#include "spaces/ublas_space.h"

//strategies
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"




//convergence criterias
#include "custom_strategies/strategies/residualbased_fluid_strategy.h"
#include "custom_strategies/strategies/residualbased_ND_fluid_strategy.h"
#include "custom_strategies/strategies/residualbased_fluid_strategy_coupled.h"
#include "custom_strategies/strategies/residualbased_lagrangian_monolithic_scheme.h"
#include "custom_strategies/strategies/residualbased_predictorcorrector_velocity_bossak_scheme.h"
#include "custom_strategies/strategies/residualbased_predictorcorrector_velocity_bossak_scheme_compressible.h"
#include "custom_strategies/strategies/residualbased_predictorcorrector_velocity_crni_scheme_compressible.h"
#include "custom_strategies/convergencecriterias/UP_criteria.h"
#include "custom_strategies/strategies/runge_kutta_fracstep_GLS_strategy.h"
#include "custom_strategies/strategies/runge_kutta_fracstep_GLS_comp_strategy.h"
#include "custom_strategies/strategies/newton_raphson_strategy.h"
//linear solvers
#include "linear_solvers/linear_solver.h"

//configuration files
#include "custom_strategies/strategies/solver_configuration.h"
#include "custom_strategies/strategies/fractionalstep_configuration.h"
#include "custom_strategies/strategies/fractional_step_strategy.h"




namespace Kratos
{

	namespace Python
	{		
		using namespace boost::python;

		void  AddCustomStrategiesToPython()
		{
			typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
			typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

			typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
			typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > BaseSolvingStrategyType;
			typedef Scheme< SparseSpaceType, LocalSpaceType > BaseSchemeType;
			typedef ConvergenceCriteria< SparseSpaceType, LocalSpaceType > TConvergenceCriteriaType;
			typedef ResidualBasedPredictorCorrectorVelocityBossakScheme< SparseSpaceType, LocalSpaceType > 					  ResidualBasedPredictorCorrectorVelocityBossakSchemeType;

			typedef ResidualBasedPredictorCorrectorVelocityBossakSchemeCompressible< SparseSpaceType, LocalSpaceType > 					  ResidualBasedPredictorCorrectorVelocityBossakSchemeCompressibleType;

			typedef ResidualBasedPredictorCorrectorVelocityCrNiSchemeCompressible< SparseSpaceType, LocalSpaceType > 					  ResidualBasedPredictorCorrectorVelocityCrNiSchemeCompressibleType;
			//********************************************************************
			//********************************************************************
			//

			class_< ResidualBasedFluidStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >,	
					bases< BaseSolvingStrategyType >,  boost::noncopyable >
				("ResidualBasedFluidStrategy", 
				init<ModelPart&, LinearSolverType::Pointer, LinearSolverType::Pointer,
				bool, bool, bool,
				double, double,
				int, int,
				unsigned int, unsigned int, unsigned int,
				bool
				>() )
				  .def("SolveStep1",&ResidualBasedFluidStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::SolveStep1)
				  .def("SolveStep2",&ResidualBasedFluidStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::SolveStep2)
				  .def("SolveStep2_Mp",&ResidualBasedFluidStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::SolveStep2_Mp)
				  .def("SolveStep3",&ResidualBasedFluidStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::SolveStep3)
				  .def("SolveStep4",&ResidualBasedFluidStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::SolveStep4)
				  .def("ActOnLonelyNodes",&ResidualBasedFluidStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::ActOnLonelyNodes)
				  .def("Clear",&ResidualBasedFluidStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::Clear)
				  .def("FractionalVelocityIteration",&ResidualBasedFluidStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::FractionalVelocityIteration)
				  .def("ConvergenceCheck",&ResidualBasedFluidStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::ConvergenceCheck)
				  .def("InitializeFractionalStep",&ResidualBasedFluidStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::InitializeFractionalStep)
				  .def("PredictVelocity",&ResidualBasedFluidStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::PredictVelocity)
				  .def("InitializeProjections",&ResidualBasedFluidStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::InitializeProjections)
				  .def("AssignInitialStepValues",&ResidualBasedFluidStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::AssignInitialStepValues)
				  .def("IterativeSolve",&ResidualBasedFluidStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::IterativeSolve)
				  .def("SavePressureIteration",&ResidualBasedFluidStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::SavePressureIteration)
				  .def("ApplyFractionalVelocityFixity",&ResidualBasedFluidStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::ApplyFractionalVelocityFixity)				;

			class_< ResidualBasedNDFluidStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >,	
					bases< BaseSolvingStrategyType >,  boost::noncopyable >
				("ResidualBasedNDFluidStrategy", 
				init<ModelPart&, LinearSolverType::Pointer, LinearSolverType::Pointer,
				bool, bool, bool,
				double, double,
				int, int,
				unsigned int, unsigned int, unsigned int,
				bool
				>() )
				  .def("SolveStep1",&ResidualBasedNDFluidStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::SolveStep1)
				  .def("SolveStep2",&ResidualBasedNDFluidStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::SolveStep2)
				  .def("SolveStep3",&ResidualBasedNDFluidStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::SolveStep3)
				  .def("SolveStep4",&ResidualBasedNDFluidStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::SolveStep4)
				  .def("ActOnLonelyNodes",&ResidualBasedNDFluidStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::ActOnLonelyNodes)
				  .def("Clear",&ResidualBasedNDFluidStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::Clear)
				  .def("FractionalVelocityIteration",&ResidualBasedNDFluidStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::FractionalVelocityIteration)
				  .def("ConvergenceCheck",&ResidualBasedNDFluidStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::ConvergenceCheck)
				  .def("InitializeFractionalStep",&ResidualBasedNDFluidStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::InitializeFractionalStep)
				  .def("PredictVelocity",&ResidualBasedNDFluidStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::PredictVelocity)
				  .def("InitializeProjections",&ResidualBasedNDFluidStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::InitializeProjections)
				  .def("AssignInitialStepValues",&ResidualBasedNDFluidStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::AssignInitialStepValues)
				  .def("IterativeSolve",&ResidualBasedNDFluidStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::IterativeSolve)
				  .def("SavePressureIteration",&ResidualBasedNDFluidStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::SavePressureIteration)
				  .def("ApplyFractionalVelocityFixity",&ResidualBasedNDFluidStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::ApplyFractionalVelocityFixity)				;	


			class_< ResidualBasedFluidStrategyCoupled< SparseSpaceType, LocalSpaceType, LinearSolverType >,	
					bases< BaseSolvingStrategyType >,  boost::noncopyable >
				("ResidualBasedFluidStrategyCoupled", 
				init<ModelPart&, LinearSolverType::Pointer, LinearSolverType::Pointer,
				bool, bool, bool,
				double, double,
				int, int,
				unsigned int, unsigned int, unsigned int,
				bool
				>() )
				  .def("SolveStep1",&ResidualBasedFluidStrategyCoupled< SparseSpaceType, LocalSpaceType, LinearSolverType >::SolveStep1)
				  .def("SolveStep2",&ResidualBasedFluidStrategyCoupled< SparseSpaceType, LocalSpaceType, LinearSolverType >::SolveStep2)
				  .def("SolveStep2_Mp",&ResidualBasedFluidStrategyCoupled< SparseSpaceType, LocalSpaceType, LinearSolverType >::SolveStep2_Mp)
				  .def("SolveStep3",&ResidualBasedFluidStrategyCoupled< SparseSpaceType, LocalSpaceType, LinearSolverType >::SolveStep3)
				  .def("SolveStep4",&ResidualBasedFluidStrategyCoupled< SparseSpaceType, LocalSpaceType, LinearSolverType >::SolveStep4)
				  .def("ActOnLonelyNodes",&ResidualBasedFluidStrategyCoupled< SparseSpaceType, LocalSpaceType, LinearSolverType >::ActOnLonelyNodes)
				  .def("Clear",&ResidualBasedFluidStrategyCoupled< SparseSpaceType, LocalSpaceType, LinearSolverType >::Clear)
				  .def("FractionalVelocityIteration",&ResidualBasedFluidStrategyCoupled< SparseSpaceType, LocalSpaceType, LinearSolverType >::FractionalVelocityIteration)
				  .def("ConvergenceCheck",&ResidualBasedFluidStrategyCoupled< SparseSpaceType, LocalSpaceType, LinearSolverType >::ConvergenceCheck)
				  .def("InitializeFractionalStep",&ResidualBasedFluidStrategyCoupled< SparseSpaceType, LocalSpaceType, LinearSolverType >::InitializeFractionalStep)
				  .def("PredictVelocity",&ResidualBasedFluidStrategyCoupled< SparseSpaceType, LocalSpaceType, LinearSolverType >::PredictVelocity)
				  .def("InitializeProjections",&ResidualBasedFluidStrategyCoupled< SparseSpaceType, LocalSpaceType, LinearSolverType >::InitializeProjections)
				  .def("AssignInitialStepValues",&ResidualBasedFluidStrategyCoupled< SparseSpaceType, LocalSpaceType, LinearSolverType >::AssignInitialStepValues)
				  .def("IterativeSolve",&ResidualBasedFluidStrategyCoupled< SparseSpaceType, LocalSpaceType, LinearSolverType >::IterativeSolve)
				  .def("SavePressureIteration",&ResidualBasedFluidStrategyCoupled< SparseSpaceType, LocalSpaceType, LinearSolverType >::SavePressureIteration)
				  .def("ApplyFractionalVelocityFixity",&ResidualBasedFluidStrategyCoupled< SparseSpaceType, LocalSpaceType, LinearSolverType >::ApplyFractionalVelocityFixity)				;	


				class_< ConvergenceCriteria< SparseSpaceType, LocalSpaceType >, boost::noncopyable >("ConvergenceCriteria", init<>() )
				.def("SetActualizeRHSFlag", &ConvergenceCriteria<SparseSpaceType, LocalSpaceType >::SetActualizeRHSFlag )
			 .def("GetActualizeRHSflag", &ConvergenceCriteria<SparseSpaceType, LocalSpaceType >::GetActualizeRHSflag )
			 .def("PreCriteria", &ConvergenceCriteria<SparseSpaceType, LocalSpaceType >::PreCriteria )
			 .def("PostCriteria", &ConvergenceCriteria<SparseSpaceType, LocalSpaceType >::PostCriteria )
			 .def("Initialize", &ConvergenceCriteria<SparseSpaceType, LocalSpaceType >::Initialize )
			 .def("InitializeSolutionStep", &ConvergenceCriteria<SparseSpaceType, LocalSpaceType >::InitializeSolutionStep )
			 .def("FinalizeSolutionStep", &ConvergenceCriteria<SparseSpaceType, LocalSpaceType >::FinalizeSolutionStep )
			 ;                  
			
			class_< UPCriteria<SparseSpaceType, LocalSpaceType >,
			         bases<ConvergenceCriteria< SparseSpaceType, LocalSpaceType > >,  
			         boost::noncopyable >
			        ("UPCriteria", init< double, double, double, double>() );


        		   class_< ResidualBasedLagrangianMonolithicScheme<SparseSpaceType,LocalSpaceType>,
        			   bases< ResidualBasedIncrementalUpdateStaticScheme<SparseSpaceType,LocalSpaceType> >,  boost::noncopyable >
             			      (
               				     "ResidualBasedLagrangianMonolithicScheme", init< int >()
                		       );

			class_< ResidualBasedPredictorCorrectorVelocityBossakSchemeType,
				bases< BaseSchemeType >,  boost::noncopyable >
					(
					"ResidualBasedPredictorCorrectorVelocityBossakScheme", init< double, double >()
					);


			class_< ResidualBasedPredictorCorrectorVelocityBossakSchemeCompressibleType,
				bases< BaseSchemeType >,  boost::noncopyable >
					(
					"ResidualBasedPredictorCorrectorVelocityBossakSchemeCompressible", init< double, double >()
					);


			class_< ResidualBasedPredictorCorrectorVelocityCrNiSchemeCompressibleType,
				bases< BaseSchemeType >,  boost::noncopyable >
					(
					"ResidualBasedPredictorCorrectorVelocityCrNiSchemeCompressible", init< double, double >()
					);

			class_< NewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >,bases< BaseSolvingStrategyType >,  boost::noncopyable >
				("NewtonRaphsonStrategy", 
				init<ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, TConvergenceCriteriaType::Pointer, int, bool, bool, bool
				>() )
				;

			class_< RungeKuttaFracStepStrategy< 2, SparseSpaceType, LocalSpaceType, LinearSolverType >,				
			
					bases< BaseSolvingStrategyType >,  boost::noncopyable >
				("RungeKuttaFracStepStrategy2D", 
				init<ModelPart&, LinearSolverType::Pointer,
				bool, bool, bool
				>() )
				  .def("SolveStep1",&RungeKuttaFracStepStrategy< 2, SparseSpaceType, LocalSpaceType, LinearSolverType >::SolveStep1)
				  .def("SolveStep2",&RungeKuttaFracStepStrategy< 2, SparseSpaceType, LocalSpaceType, LinearSolverType >::SolveStep2)
				  .def("SolveStep3",&RungeKuttaFracStepStrategy< 2, SparseSpaceType, LocalSpaceType, LinearSolverType >::SolveStep3)
				  
				  .def("Clear",&ResidualBasedFluidStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::Clear);

			class_< RungeKuttaFracStepStrategy< 3, SparseSpaceType, LocalSpaceType, LinearSolverType >,	
					bases< BaseSolvingStrategyType >,  boost::noncopyable >
				("RungeKuttaFracStepStrategy3D", 
				init<ModelPart&, LinearSolverType::Pointer,
				bool, bool, bool
				>() )
				  .def("SolveStep1",&RungeKuttaFracStepStrategy< 3, SparseSpaceType, LocalSpaceType, LinearSolverType >::SolveStep1)
				  .def("SolveStep2",&RungeKuttaFracStepStrategy< 3, SparseSpaceType, LocalSpaceType, LinearSolverType >::SolveStep2)
				  .def("SolveStep3",&RungeKuttaFracStepStrategy< 3, SparseSpaceType, LocalSpaceType, LinearSolverType >::SolveStep3)
				  
				  .def("Clear",&ResidualBasedFluidStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::Clear);

			
			class_< RungeKuttaFracStepCompStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >,	
					bases< BaseSolvingStrategyType >,  boost::noncopyable >
				("RungeKuttaFracStepCompStrategy", 
				init<ModelPart&, LinearSolverType::Pointer,
				bool, bool, bool,
				int 
				>() )
				  .def("SolveStep1",&RungeKuttaFracStepCompStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::SolveStep1)
				  .def("SolveStep2",&RungeKuttaFracStepCompStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::SolveStep2)
				  .def("SolveStep3",&RungeKuttaFracStepCompStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::SolveStep3)
				  
				  .def("Clear",&ResidualBasedFluidStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::Clear);



                        //********************************************************************************************
                       class_< SolverConfiguration<SparseSpaceType, LocalSpaceType, LinearSolverType >,
			         boost::noncopyable >
			        ("SolverConfiguration", init< ModelPart&, unsigned int>() )
                           .def("GetActualizeRHSflag", &ConvergenceCriteria<SparseSpaceType, LocalSpaceType >::GetActualizeRHSflag )
                           .def("GetActualizeRHSflag", &ConvergenceCriteria<SparseSpaceType, LocalSpaceType >::GetActualizeRHSflag )
                           ;

                       class_< FractionalStepConfiguration<SparseSpaceType, LocalSpaceType, LinearSolverType >,
                               bases< SolverConfiguration<SparseSpaceType, LocalSpaceType, LinearSolverType > >,
			         boost::noncopyable >
			        ("FractionalStepConfiguration", init< ModelPart&, LinearSolverType::Pointer, LinearSolverType::Pointer, 
                                                                        unsigned int, unsigned int, bool >() );



                        //********************************************************************************************
			class_< FractionalStepStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >,
					bases< BaseSolvingStrategyType >,  boost::noncopyable >
				("FractionalStepStrategy",
				init<ModelPart&,
                                    SolverConfiguration<SparseSpaceType, LocalSpaceType, LinearSolverType >&,
                                    bool,
                                    double, double,
                                    int, int,
                                    unsigned int, unsigned int,
                                    bool
                                    >() )
				  .def("SolveStep1",&FractionalStepStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::SolveStep1)
				  .def("SolveStep2",&FractionalStepStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::SolveStep2)
				  .def("SolveStep3",&FractionalStepStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::SolveStep3)
				  .def("SolveStep4",&FractionalStepStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::SolveStep4)
				  .def("ActOnLonelyNodes",&FractionalStepStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::ActOnLonelyNodes)
				  .def("Clear",&FractionalStepStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::Clear)
				  .def("FractionalVelocityIteration",&FractionalStepStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::FractionalVelocityIteration)
				  .def("ConvergenceCheck",&FractionalStepStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::ConvergenceCheck)
				  .def("InitializeFractionalStep",&FractionalStepStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::InitializeFractionalStep)
				  .def("PredictVelocity",&FractionalStepStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::PredictVelocity)
				  .def("InitializeProjections",&FractionalStepStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::InitializeProjections)
				  .def("AssignInitialStepValues",&FractionalStepStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::AssignInitialStepValues)
				  .def("IterativeSolve",&FractionalStepStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::IterativeSolve)
				  .def("SavePressureIteration",&FractionalStepStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::SavePressureIteration)
				  .def("ApplyFractionalVelocityFixity",&FractionalStepStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::ApplyFractionalVelocityFixity)				;





		}

	}  // namespace Python.

} // Namespace Kratos

