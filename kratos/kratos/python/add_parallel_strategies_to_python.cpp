/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

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
//   Last modified by:    $Author: nelson $
//   Date:                $Date: 2008-12-09 15:23:36 $
//   Revision:            $Revision: 1.6 $
//
//


// System includes 

//nothing will be compiled if an openmp compiler is not found

// External includes 
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/timer.hpp> 
#include "python/add_parallel_strategies_to_python.h"
#include <cstring>

#ifdef _OPENMP

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "spaces/parallel_ublas_space.h"
#include "spaces/ublas_space.h"

//strategies
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"

//schemes
#include "solving_strategies/schemes/scheme.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"

//convergence criterias
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "solving_strategies/convergencecriterias/displacement_criteria.h"
//#include "solving_strategies/convergencecriterias/residual_criteria.h"
//#include "solving_strategies/convergencecriterias/and_criteria.h"

//Builder And Solver
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "solving_strategies/builder_and_solvers/parallel_residualbased_elimination_builder_and_solver.h"
#include "solving_strategies/builder_and_solvers/parallel_residualbased_elimination_builder_and_solver_deactivation.h"


//linear solvers
#include "linear_solvers/linear_solver.h"

//utilities
#include "python/pointer_vector_set_python_interface.h"

#endif


namespace Kratos
{

	namespace Python
	{		

		using namespace boost::python;


#ifdef _OPENMP //nothing will be compiled if an openmp compiler is not found


		typedef ParallelUblasSpace<double, CompressedMatrix, Vector> ParallelSparseSpaceType;
		typedef UblasSpace<double, Matrix, Vector> ParallelLocalSpaceType;

		void ParallelResizeMatrix(ParallelSparseSpaceType& dummy, ParallelSparseSpaceType::MatrixType& A, unsigned int i1, unsigned int i2)
		{	dummy.Resize(A,i1,i2);	}			

		void ParallelResizeVector(ParallelSparseSpaceType& dummy, ParallelSparseSpaceType::VectorType& x, unsigned int i1)
		{	dummy.Resize(x,i1);	}

		void ParallelSetToZeroMatrix(ParallelSparseSpaceType& dummy, ParallelSparseSpaceType::MatrixType& A)
		{	dummy.SetToZero(A);	}

		void ParallelSetToZeroVector(ParallelSparseSpaceType& dummy, ParallelSparseSpaceType::VectorType& x)
		{	dummy.SetToZero(x);	}

		void ParallelClearMatrix(ParallelSparseSpaceType& dummy, ParallelSparseSpaceType::MatrixPointerType& A)
		{	dummy.Clear(A);	}

		void ParallelClearVector(ParallelSparseSpaceType& dummy, ParallelSparseSpaceType::VectorPointerType& x)
		{	dummy.Clear(x);	}

		double ParallelTwoNorm(ParallelSparseSpaceType& dummy, ParallelSparseSpaceType::VectorType& x)
		{	return dummy.TwoNorm(x);	}

		void ParallelMoveMesh( Scheme< ParallelSparseSpaceType, ParallelLocalSpaceType >& dummy, ModelPart::NodesContainerType& rNodes)
		{	
				for(ModelPart::NodeIterator i = rNodes.begin() ; i != rNodes.end() ; ++i)
				{
					const array_1d<double,3>& disp = i->FastGetSolutionStepValue(DISPLACEMENT);
					(i)->X() = (i)->X0() + disp[0];
					(i)->Y() = (i)->Y0() + disp[1];
					(i)->Z() = (i)->Z0() + disp[2];
				}
		}

		ParallelSparseSpaceType::MatrixPointerType CreateEmptyMatrixPointer(ParallelSparseSpaceType& dummy){ return dummy.CreateEmptyMatrixPointer();}
		ParallelSparseSpaceType::VectorPointerType CreateEmptyVectorPointer(ParallelSparseSpaceType& dummy){ return dummy.CreateEmptyVectorPointer();}
	
#endif //nothing will be compiled if an openmp compiler is not found


		void  AddParallelStrategiesToPython()
		{
//nothing will be compiled if an openmp compiler is not found
#ifdef _OPENMP
//			typedef ParallelUblasSpace<double, CompressedMatrix, Vector> ParallelSparseSpaceType;
//			typedef UblasSpace<double, Matrix, Vector> ParallelLocalSpaceType;

			typedef LinearSolver<ParallelSparseSpaceType, ParallelLocalSpaceType > ParallelLinearSolverType;
			typedef SolvingStrategy< ParallelSparseSpaceType, ParallelLocalSpaceType, ParallelLinearSolverType > ParallelBaseSolvingStrategyType;
			typedef Scheme< ParallelSparseSpaceType, ParallelLocalSpaceType > ParallelBaseSchemeType;

			//********************************************************************
			//********************************************************************
			//strategy base class
			class_< ParallelBaseSolvingStrategyType, boost::noncopyable >("ParallelSolvingStrategy", init< ModelPart&, bool >() )
				.def("Predict", &ParallelBaseSolvingStrategyType::Predict )
				.def("Solve", &ParallelBaseSolvingStrategyType::Solve )
				.def("IsConverged", &ParallelBaseSolvingStrategyType::IsConverged )
				.def("CalculateOutputData", &ParallelBaseSolvingStrategyType::CalculateOutputData )
				.def("SetEchoLevel", &ParallelBaseSolvingStrategyType::SetEchoLevel )
				.def("GetEchoLevel", &ParallelBaseSolvingStrategyType::GetEchoLevel )
				.def("SetRebuildLevel", &ParallelBaseSolvingStrategyType::SetRebuildLevel )
				.def("GetRebuildLevel", &ParallelBaseSolvingStrategyType::GetRebuildLevel )
				.def("SetMoveMeshFlag", &ParallelBaseSolvingStrategyType::SetMoveMeshFlag )
				.def("MoveMeshFlag", &ParallelBaseSolvingStrategyType::MoveMeshFlag )
				.def("MoveMesh", &ParallelBaseSolvingStrategyType::MoveMesh )
				.def("Clear", &ParallelBaseSolvingStrategyType::Clear )
				//.def("GetModelPart", &ParallelBaseSolvingStrategyType::GetModelPart )
				; 

			
			class_< ResidualBasedLinearStrategy< ParallelSparseSpaceType, ParallelLocalSpaceType, ParallelLinearSolverType >,bases< ParallelBaseSolvingStrategyType >,  boost::noncopyable >
				("ParallelResidualBasedLinearStrategy", 
				init<ModelPart&,ParallelBaseSchemeType::Pointer, ParallelLinearSolverType::Pointer, bool, bool, bool, bool	>() )
				;

			typedef ConvergenceCriteria< ParallelSparseSpaceType, ParallelLocalSpaceType > TConvergenceCriteriaType;

			class_< ResidualBasedNewtonRaphsonStrategy< ParallelSparseSpaceType, ParallelLocalSpaceType, ParallelLinearSolverType >,bases< ParallelBaseSolvingStrategyType >,  boost::noncopyable >
				("ParallelResidualBasedNewtonRaphsonStrategy", 
				init<ModelPart&, ParallelBaseSchemeType::Pointer, ParallelLinearSolverType::Pointer, TConvergenceCriteriaType::Pointer, int, bool, bool, bool
				>() )
				;


			//********************************************************************
			//********************************************************************
			class_< ParallelBaseSchemeType, boost::noncopyable >
			("ParallelScheme", init< >() )
			.def("Initialize", &ParallelBaseSchemeType::Initialize )
			.def("SchemeIsInitialized", &ParallelBaseSchemeType::SchemeIsInitialized )
			.def("ElementsAreInitialized", &ParallelBaseSchemeType::ElementsAreInitialized )
			.def("InitializeElements", &ParallelBaseSchemeType::InitializeElements )
			.def("InitializeSolutionStep", &ParallelBaseSchemeType::InitializeSolutionStep )
			.def("FinalizeSolutionStep", &ParallelBaseSchemeType::FinalizeSolutionStep )
			.def("InitializeNonLinIteration", &ParallelBaseSchemeType::InitializeNonLinIteration )
			.def("FinalizeNonLinIteration", &ParallelBaseSchemeType::FinalizeNonLinIteration )
			.def("Predict", &ParallelBaseSchemeType::Predict )
			.def("Update", &ParallelBaseSchemeType::Update )
			.def("CalculateOutputData", &ParallelBaseSchemeType::CalculateOutputData )
			.def("Clean", &ParallelBaseSchemeType::Clean )
			.def("MoveMesh", ParallelMoveMesh )
			;

			class_< ResidualBasedIncrementalUpdateStaticScheme< ParallelSparseSpaceType, ParallelLocalSpaceType>,	
					bases< ParallelBaseSchemeType >,  boost::noncopyable >
				(
					"ParallelResidualBasedIncrementalUpdateStaticScheme", init< >() 
				);


			//********************************************************************
			//********************************************************************
			//convergence criteria base class
			 class_< ConvergenceCriteria< ParallelSparseSpaceType, ParallelLocalSpaceType >, boost::noncopyable >("ParallelConvergenceCriteria", init<>() )
				.def("SetActualizeRHSFlag", &ConvergenceCriteria<ParallelSparseSpaceType, ParallelLocalSpaceType >::SetActualizeRHSFlag )
			 .def("GetActualizeRHSflag", &ConvergenceCriteria<ParallelSparseSpaceType, ParallelLocalSpaceType >::GetActualizeRHSflag )
			 .def("PreCriteria", &ConvergenceCriteria<ParallelSparseSpaceType, ParallelLocalSpaceType >::PreCriteria )
			 .def("PostCriteria", &ConvergenceCriteria<ParallelSparseSpaceType, ParallelLocalSpaceType >::PostCriteria )
			 .def("Initialize", &ConvergenceCriteria<ParallelSparseSpaceType, ParallelLocalSpaceType >::Initialize )
			 .def("InitializeSolutionStep", &ConvergenceCriteria<ParallelSparseSpaceType, ParallelLocalSpaceType >::InitializeSolutionStep )
			 .def("FinalizeSolutionStep", &ConvergenceCriteria<ParallelSparseSpaceType, ParallelLocalSpaceType >::FinalizeSolutionStep )
			 ;                  
			
			class_< DisplacementCriteria<ParallelSparseSpaceType, ParallelLocalSpaceType >,
			         bases<ConvergenceCriteria< ParallelSparseSpaceType, ParallelLocalSpaceType > >,  
			         boost::noncopyable >
			        ("ParallelDisplacementCriteria", init< double, double>() );
			
/*			class_< ResidualCriteria< ParallelSparseSpaceType >,
			         bases<ConvergenceCriteria< ParallelSparseSpaceType > >,  
			         boost::noncopyable >
			        ("ResidualCriteria", init<Model::Pointer, double >() );
			
			class_< AndCriteria< ParallelSparseSpaceType >,
			         bases<ConvergenceCriteria< ParallelSparseSpaceType > >,  
			         boost::noncopyable >
			        ("AndCriteria", init<Model::Pointer, ConvergenceCriteria< ParallelSparseSpaceType >::Pointer, ConvergenceCriteria< ParallelSparseSpaceType >::Pointer >()*/
			//);

			//********************************************************************
			//********************************************************************
			//
			
			//Builder and Solver
			typedef BuilderAndSolver< ParallelSparseSpaceType, ParallelLocalSpaceType, ParallelLinearSolverType > ParallelBuilderAndSolverType;

			class_< ParallelBuilderAndSolverType::DofsArrayType, boost::noncopyable >("ParallelDofsArrayType",	init<>() );
			 
			class_< ParallelBuilderAndSolverType, boost::noncopyable >("ParallelBuilderAndSolver",	init<ParallelLinearSolverType::Pointer>() )
			 .def("SetCalculateReactionsFlag", &ParallelBuilderAndSolverType::SetCalculateReactionsFlag )
			 .def("GetCalculateReactionsFlag", &ParallelBuilderAndSolverType::GetCalculateReactionsFlag )
			 .def("SetDofSetIsInitializedFlag", &ParallelBuilderAndSolverType::SetDofSetIsInitializedFlag )
			 .def("GetDofSetIsInitializedFlag", &ParallelBuilderAndSolverType::GetDofSetIsInitializedFlag )
			 .def("SetReshapeMatrixFlag", &ParallelBuilderAndSolverType::SetReshapeMatrixFlag )
			 .def("GetReshapeMatrixFlag", &ParallelBuilderAndSolverType::GetReshapeMatrixFlag )
			 .def("GetEquationSystemSize", &ParallelBuilderAndSolverType::GetEquationSystemSize )
			 .def("BuildLHS", &ParallelBuilderAndSolverType::BuildLHS )
			 .def("BuildRHS", &ParallelBuilderAndSolverType::BuildRHS )
			 .def("Build", &ParallelBuilderAndSolverType::Build )
			 .def("SystemSolve", &ParallelBuilderAndSolverType::SystemSolve )
			 .def("BuildAndSolve", &ParallelBuilderAndSolverType::BuildAndSolve )
			 .def("BuildRHSAndSolve", &ParallelBuilderAndSolverType::BuildRHSAndSolve )
			 .def("ApplyDirichletConditions", &ParallelBuilderAndSolverType::ApplyDirichletConditions )
			 .def("SetUpDofSet", &ParallelBuilderAndSolverType::SetUpDofSet )
			 .def("GetDofSet", &ParallelBuilderAndSolverType::GetDofSet, return_internal_reference<>() )
			 .def("SetUpSystem", &ParallelBuilderAndSolverType::SetUpSystem )
			 .def("ResizeAndInitializeVectors", &ParallelBuilderAndSolverType::ResizeAndInitializeVectors )
			 .def("InitializeSolutionStep", &ParallelBuilderAndSolverType::InitializeSolutionStep )
			 .def("FinalizeSolutionStep", &ParallelBuilderAndSolverType::FinalizeSolutionStep )
			 .def("CalculateReactions", &ParallelBuilderAndSolverType::CalculateReactions )
			 .def("Clear", &ParallelBuilderAndSolverType::Clear )
			 .def("SetEchoLevel", &ParallelBuilderAndSolverType::SetEchoLevel )
			 .def("GetEchoLevel", &ParallelBuilderAndSolverType::GetEchoLevel )
			 ;

			typedef ParallelResidualBasedEliminationBuilderAndSolver< ParallelSparseSpaceType, ParallelLocalSpaceType, ParallelLinearSolverType > ParallelResidualBasedEliminationBuilderAndSolverType;

			class_< ParallelResidualBasedEliminationBuilderAndSolverType, bases<ParallelBuilderAndSolverType>, boost::noncopyable> ("ParallelResidualBasedEliminationBuilderAndSolver", init< ParallelLinearSolverType::Pointer>() );
            
            typedef ParallelResidualBasedEliminationBuilderAndSolverDeactivation< ParallelSparseSpaceType, ParallelLocalSpaceType, ParallelLinearSolverType > ParallelResidualBasedEliminationBuilderAndSolverDeactivationType;
            
            class_< ParallelResidualBasedEliminationBuilderAndSolverDeactivationType, bases<ParallelBuilderAndSolverType>, boost::noncopyable> ("ParallelResidualBasedEliminationBuilderAndSolverDeactivation", init< ParallelLinearSolverType::Pointer>() );


			//********************************************************************
			//********************************************************************

			class_< ParallelSparseSpaceType, boost::noncopyable >("ParallelUblasSparseSpace",	init<>() )
			 .def("ClearMatrix", ParallelClearMatrix )
			 .def("ClearVector", ParallelClearVector )
			 .def("ResizeMatrix", ParallelResizeMatrix )
			 .def("ResizeVector", ParallelResizeVector )
			 .def("SetToZeroMatrix", ParallelSetToZeroMatrix )
			 .def("SetToZeroVector", ParallelSetToZeroVector )
			 .def("TwoNorm", ParallelTwoNorm )
			 .def("CreateEmptyMatrixPointer",CreateEmptyMatrixPointer)
			 .def("CreateEmptyVectorPointer",CreateEmptyVectorPointer)
			;
#endif
}


	}  // namespace Python.

} // Namespace Kratos

