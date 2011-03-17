/*
==============================================================================
KratosTestApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel 
pooyan@cimne.upc.edu 
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


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
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/timer.hpp> 


// Project includes
#include "includes/define.h"
#include "custom_python/add_custom_strategies_to_python.h"

#include "spaces/ublas_space.h"

//strategies
#include "solving_strategies/strategies/solving_strategy.h"
#include "custom_strategies/residualbased_elimination_builder_and_solver_pod.h"
#include "custom_strategies/residualbased_elimination_builder_and_solver_pod_WithPressure.h"
#include "custom_strategies/residualbased_elimination_builder_and_solver_standard.h"
#include "custom_strategies/residualbased_elimination_builder_and_solver_standard_biggerlocalM.h"
#include "custom_strategies/residualbased_elimination_builder_and_solver_standard_biggerlocalM_NoP.h"

//linear solvers
#include "linear_solvers/linear_solver.h"



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
                    
			typedef BuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType>
				BuilderAndSolverType;
			//********************************************************************
			//********************************************************************
// 			class_< TestStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >,	
// 					bases< BaseSolvingStrategyType >,  boost::noncopyable >
// 				("TestStrategy", 
// 				init<ModelPart&, LinearSolverType::Pointer, int, int, bool >() )
// 				.def("MoveNodes",&TestStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::MoveNodes)
// 				;
			typedef ResidualBasedEliminationBuilderAndSolverStandard< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedEliminationBuilderAndSolverStandardType;
			class_< ResidualBasedEliminationBuilderAndSolverStandardType, bases<BuilderAndSolverType>, boost::noncopyable> ("ResidualBasedEliminationBuilderAndSolverStandard", init< LinearSolverType::Pointer>() );

			typedef ResidualBasedEliminationBuilderAndSolverStandard_biggerlocalM< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedEliminationBuilderAndSolverBiggerPODType;
			class_< ResidualBasedEliminationBuilderAndSolverBiggerPODType, bases<BuilderAndSolverType>, boost::noncopyable> ("ResidualBasedEliminationBuilderAndSolverStandard_BiggerLocalM", init< LinearSolverType::Pointer>() );
			
			typedef ResidualBasedEliminationBuilderAndSolverPOD< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedEliminationBuilderAndSolverPODType;
			class_< ResidualBasedEliminationBuilderAndSolverPODType, bases<BuilderAndSolverType>, boost::noncopyable> ("ResidualBasedEliminationBuilderAndSolverPOD", init< LinearSolverType::Pointer>() );
			
			typedef ResidualBasedEliminationBuilderAndSolverPOD_WithPressure< SparseSpaceType, LocalSpaceType, LinearSolverType >    ResidualBasedEliminationBuilderAndSolverPOD_WithPressureType;
			class_< ResidualBasedEliminationBuilderAndSolverPOD_WithPressureType, bases<BuilderAndSolverType>, boost::noncopyable> ("ResidualBasedEliminationBuilderAndSolverPOD_WithPressure", init< LinearSolverType::Pointer>() );
																		 
			typedef ResidualBasedEliminationBuilderAndSolverStandard_biggerlocalM_NoP< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedEliminationBuilderAndSolverStandard_biggerlocalM_NoPType;
			class_< ResidualBasedEliminationBuilderAndSolverStandard_biggerlocalM_NoPType, bases<BuilderAndSolverType>, boost::noncopyable> ("ResidualBasedEliminationBuilderAndSolverStandard_biggerlocalM_NoP", init< LinearSolverType::Pointer>() );
		}

	}  // namespace Python.

} // Namespace Kratos

