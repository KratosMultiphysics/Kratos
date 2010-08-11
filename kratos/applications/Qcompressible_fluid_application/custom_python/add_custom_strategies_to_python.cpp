/*
==============================================================================
KratosR1IncompressibleFluidApplication 
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
//   Last modified by:    $Author: jmarti $
//   Date:                $Date: 2009-01-23 14:33:59 $
//   Revision:            $Revision: 1.1 $
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
#include "custom_strategies/strategies/residualbasedq_fluid_strategy.h"
#include "custom_strategies/strategies/residualbased_predictorcorrector_velocity_bossak_scheme.h"


//linear solvers
#include "linear_solvers/linear_solver.h"

//convergence criterias

#include "custom_strategies/convergencecriterias/UP_criteria.h"

//linear solvers

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
			typedef ResidualBasedPredictorCorrectorVelocityBossakScheme< SparseSpaceType, LocalSpaceType > 					  ResidualBasedPredictorCorrectorVelocityBossakSchemeType;
			//********************************************************************
			//********************************************************************
			//

			
			class_< UPCriteria<SparseSpaceType, LocalSpaceType >, bases<ConvergenceCriteria< SparseSpaceType, LocalSpaceType > >,  boost::noncopyable >
			("UPCriteria", init< double, double, double, double>() );
			
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
	  
	  class_< ResidualBasedPredictorCorrectorVelocityBossakSchemeType,
		  bases< BaseSchemeType >,  boost::noncopyable >
	  (
	   "ResidualBasedPredictorCorrectorVelocityBossakScheme", init< double, double >()
					);

			


			/*class_< ResidualBasedNDFluidStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >,	
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
				  .def("ApplyFractionalVelocityFixity",&ResidualBasedNDFluidStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::ApplyFractionalVelocityFixity)				;	*/


			/*class_< ResidualBasedFluidStrategyCoupled< SparseSpaceType, LocalSpaceType, LinearSolverType >,	
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
				  .def("ApplyFractionalVelocityFixity",&ResidualBasedFluidStrategyCoupled< SparseSpaceType, LocalSpaceType, LinearSolverType >::ApplyFractionalVelocityFixity)				;*/	














		}

	}  // namespace Python.

} // Namespace Kratos

