//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//
// System includes


// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define.h"
#include "custom_python/add_custom_strategies_to_python.h"

#include "spaces/ublas_space.h"

//strategies
#include "solving_strategies/strategies/solving_strategy.h"
// #include "custom_strategies/hexahedra_newton_raphson_strategy.h"
#include "custom_strategies/femdem_residual_criteria.h"

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
		typedef Scheme< SparseSpaceType, LocalSpaceType > BaseSchemeType;
		typedef FemDemResidualCriteria< SparseSpaceType,  LocalSpaceType > FemDemResidualCriteriaType;

		// Base types
		typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
		typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > BaseSolvingStrategyType;
		typedef ConvergenceCriteria< SparseSpaceType, LocalSpaceType > ConvergenceCriteriaType;
		// typedef HexahedraNewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > HexahedraNewtonRaphsonStrategyType;

		// class_< HexahedraNewtonRaphsonStrategyType,
		// 		typename HexahedraNewtonRaphsonStrategyType::Pointer,
		// 		BaseSolvingStrategyType  >  (m, "HexahedraNewtonRaphsonStrategy")
		// 		.def(init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, int, bool, bool, bool >())
		// 		;
	
		class_<FemDemResidualCriteria<SparseSpaceType, LocalSpaceType >,
			typename FemDemResidualCriteria<SparseSpaceType, LocalSpaceType >::Pointer,
			ConvergenceCriteriaType >
			(m,"FemDemResidualCriteria")
			.def(init< double, double>())
			;
	}

}  // namespace Python.

} // Namespace Kratos
