//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Antonia Larese 
//


// System includes


// External includes
// #include <boost/timer.hpp>


// Project includes
#include "includes/define_python.h"
#include "custom_python/add_custom_strategies_to_python.h"

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

		void  AddCustomStrategiesToPython(pybind11::module& pymodule)
		{
			typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
			typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

			typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
			typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > BaseSolvingStrategyType;
			typedef Scheme< SparseSpaceType, LocalSpaceType > BaseSchemeType;

			//********************************************************************
			//********************************************************************
// 			class_< TestStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >,
// 					bases< BaseSolvingStrategyType >,  boost::noncopyable >
// 				("TestStrategy",
// 				init<ModelPart&, LinearSolverType::Pointer, int, int, bool >() )
// 				.def("MoveNodes",&TestStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::MoveNodes)
// 				;

		}

	}  // namespace Python.

} // Namespace Kratos
