//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Author Julio Marti
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
//#include "custom_strategies/strategies/residualbased_convdiffr_strategy.h"
#include "custom_strategies/strategies/residualbased_convdiffr_strategy_nonlinear.h"

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

			//********************************************************************
			//********************************************************************
			//


/*			class_< ResidualBasedConvectionDiffusionrStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >,	
					bases< BaseSolvingStrategyType >,  boost::noncopyable >
				("ResidualBasedConvectionDiffusionrStrategy", 
				init<	ModelPart&, LinearSolverType::Pointer,	bool, int, int	>() )
				  .def("Clear",&ResidualBasedConvectionDiffusionrStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::Clear)
				;
*/
			class_< ResidualBasedConvectionDiffusionrStrategyNonLinear< SparseSpaceType, LocalSpaceType, LinearSolverType >,	
					bases< BaseSolvingStrategyType >,  boost::noncopyable >
				("ResidualBasedConvectionDiffusionrStrategyNonLinear", 
				init<	ModelPart&, LinearSolverType::Pointer,	bool, int, int ,double	>() )
				  .def("Clear",&ResidualBasedConvectionDiffusionrStrategyNonLinear< SparseSpaceType, LocalSpaceType, LinearSolverType >::Clear)
				;	
		}

	}  // namespace Python.

} // Namespace Kratos

