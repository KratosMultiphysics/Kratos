/*
==============================================================================
KratosStructuralApplication 
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
//   Last modified by:    $Author: janosch $
//   Date:                $Date: 2008-10-23 12:48:56 $
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
#include "custom_python/add_custom_parallel_strategies_to_python.h"

#ifdef _OPENMP //nothing will be compiled if an openmp compiler is not found

	#include "spaces/ublas_space.h"
	#include "spaces/parallel_ublas_space.h"

	//convergence criteria
	#include "solving_strategies/convergencecriterias/convergence_criteria.h"

	//strategies
	#include "solving_strategies/strategies/solving_strategy.h"

	//schemes
	#include "solving_strategies/schemes/scheme.h"
	#include "custom_strategies/schemes/residualbased_predictorcorrector_bossak_scheme.h"
	#include "custom_strategies/schemes/residualbased_newmark_scheme.h"

	//builder_and_solvers
	#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
	#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"

	//linear solvers
	#include "linear_solvers/linear_solver.h"
#endif

namespace Kratos
{

	namespace Python
	{		
		using namespace boost::python;

		void  AddCustomParallelStrategiesToPython()
		{
#ifdef _OPENMP //nothing will be compiled if an openmp compiler is not found
			typedef UblasSpace<double, Matrix, Vector> ParallelLocalSpaceType;
			typedef ParallelUblasSpace<double, CompressedMatrix, Vector> ParallelSparseSpaceType;
//			typedef UblasSpace<double, CompressedMatrix, Vector> ParallelSparseSpaceType;
//std::cout << "ATTENTION Ublas space used instead of Parallel Ublas Space" << std::endl;

			typedef LinearSolver<ParallelSparseSpaceType, ParallelLocalSpaceType > ParallelLinearSolverType;
			typedef SolvingStrategy< ParallelSparseSpaceType, ParallelLocalSpaceType, ParallelLinearSolverType > ParallelBaseSolvingStrategyType;
			typedef Scheme< ParallelSparseSpaceType, ParallelLocalSpaceType > ParallelBaseSchemeType;

			

			typedef ResidualBasedPredictorCorrectorBossakScheme< ParallelSparseSpaceType, ParallelLocalSpaceType > ParallelResidualBasedPredictorCorrectorBossakSchemeType;
            typedef ResidualBasedNewmarkScheme<ParallelSparseSpaceType, ParallelLocalSpaceType> ParallelResidualBasedNewmarkSchemeType;
					
			typedef ConvergenceCriteria< ParallelSparseSpaceType, ParallelLocalSpaceType > ParallelConvergenceCriteriaBaseType;




			//********************************************************************
            		//********************************************************************

			class_< ParallelResidualBasedPredictorCorrectorBossakSchemeType,
			bases< ParallelBaseSchemeType >,  boost::noncopyable >
					(
					"ParallelResidualBasedPredictorCorrectorBossakScheme", init< double >()
					);
            
            class_< ParallelResidualBasedNewmarkSchemeType,
            bases<ParallelBaseSchemeType >, boost::noncopyable >
                    (
                     "ParallelResidualBasedNewmarkScheme",init<double>()
                    );
#endif
			
		}

	}  // namespace Python.

} // Namespace Kratos

