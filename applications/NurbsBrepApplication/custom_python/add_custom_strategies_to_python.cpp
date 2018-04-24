//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: nurbs_brep_application/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes


// External includes


// Project includes
#include "includes/define_python.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "includes/define.h"

#include "spaces/ublas_space.h"

//strategies
#include "solving_strategies/strategies/solving_strategy.h"

//linear solvers
#include "linear_solvers/linear_solver.h"



namespace Kratos
{

	namespace Python
	{
		using namespace pybind11;

		void  AddCustomStrategiesToPython(pybind11::module& m)
		{
			typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
			typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

			typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
			typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > BaseSolvingStrategyType;
			typedef Scheme< SparseSpaceType, LocalSpaceType > BaseSchemeType;


		}

	}  // namespace Python.

} // Namespace Kratos
