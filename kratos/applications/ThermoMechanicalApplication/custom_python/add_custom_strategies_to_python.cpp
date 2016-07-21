//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Jordi Rubio
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
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "custom_strategies/residualbased_incrementalupdate_variable_property_static_scheme.h"
#include "custom_strategies/general_residualbased_incrementalupdate_variable_property_static_scheme.h"

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

//     typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
    //typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > BaseSolvingStrategyType;
    //typedef Scheme< SparseSpaceType, LocalSpaceType > BaseSchemeType;

    //********************************************************************
    //********************************************************************
    class_< ResidualBasedIncrementalUpdateStaticVariablePropertyScheme< SparseSpaceType, LocalSpaceType>,
	    bases< ResidualBasedIncrementalUpdateStaticScheme <SparseSpaceType, LocalSpaceType> >, boost::noncopyable >
	    (
	    "ResidualBasedIncrementalUpdateStaticVariablePropertyScheme", init< >()
	    );

    class_< GeneralResidualBasedIncrementalUpdateStaticVariablePropertyScheme< SparseSpaceType, LocalSpaceType>,
	    bases< ResidualBasedIncrementalUpdateStaticScheme <SparseSpaceType, LocalSpaceType> >, boost::noncopyable >
	    (
	    "GeneralResidualBasedIncrementalUpdateStaticVariablePropertyScheme", init< >()
	    );

// 			class_< TestStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >,
// 					bases< BaseSolvingStrategyType >,  boost::noncopyable >
// 				("TestStrategy",
// 				init<ModelPart&, LinearSolverType::Pointer, int, int, bool >() )
// 				.def("MoveNodes",&TestStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::MoveNodes)
// 				;

}

}  // namespace Python.

} // Namespace Kratos

