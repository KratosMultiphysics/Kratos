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
//#include <boost/python.hpp>
//#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
//#include <boost/timer.hpp>
#include <pybind11/pybind11.h>

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
		namespace py = pybind11;

		void  AddCustomStrategiesToPython(pybind11::module& m)
		{
			typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
			typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

			typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
			typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > BaseSolvingStrategyType;
			typedef Scheme< SparseSpaceType, LocalSpaceType > BaseSchemeType;

			//********************************************************************
			//********************************************************************
			//


		py::class_< ResidualBasedConvectionDiffusionrStrategyNonLinear< SparseSpaceType, LocalSpaceType, LinearSolverType >,ResidualBasedConvectionDiffusionrStrategyNonLinear< SparseSpaceType, LocalSpaceType, LinearSolverType >::Pointer,BaseSolvingStrategyType> (m,"ResidualBasedConvectionDiffusionrStrategyNonLinear").def(init<ModelPart&, LinearSolverType::Pointer,	bool, int, int ,double	>() );


		}

	}  // namespace Python.

} // Namespace Kratos

