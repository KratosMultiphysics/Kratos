// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___ 
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _| 
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License: BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:  Riccardo Rossi
//

// System includes


// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define_python.h"
#include "custom_python/add_custom_strategies_to_python.h"

#include "spaces/ublas_space.h"

//strategies
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "custom_strategies/strategies/residualbased_eulerian_convdiff_strategy.h"
#include "custom_strategies/strategies/residualbased_semi_eulerian_convdiff_strategy.h"


//convergence criterias
#include "custom_strategies/strategies/residualbased_convdiff_strategy.h"
#include "custom_strategies/strategies/residualbased_convdiff_strategy_nonlinear.h"

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
//     typedef Scheme< SparseSpaceType, LocalSpaceType > BaseSchemeType;

    //********************************************************************
    //********************************************************************
    //


    class_< ResidualBasedConvectionDiffusionStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >,
            ResidualBasedConvectionDiffusionStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::Pointer,
            BaseSolvingStrategyType >
            (m,"ResidualBasedConvectionDiffusionStrategy")
            .def(init<	ModelPart&, LinearSolverType::Pointer,	bool, int, int	>() )
            .def("Clear",&ResidualBasedConvectionDiffusionStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::Clear)
            ;
    
    class_< ResidualBasedEulerianConvectionDiffusionStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >,
            ResidualBasedEulerianConvectionDiffusionStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::Pointer,
            BaseSolvingStrategyType >
            (m,"ResidualBasedEulerianConvectionDiffusionStrategy")
            .def( init<	ModelPart&, LinearSolverType::Pointer,	bool, int	>() )
            .def("Clear",&ResidualBasedEulerianConvectionDiffusionStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::Clear)
            ;
            
    class_< ResidualBasedSemiEulerianConvectionDiffusionStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >,
            ResidualBasedSemiEulerianConvectionDiffusionStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::Pointer,
            BaseSolvingStrategyType >
            (m,"ResidualBasedSemiEulerianConvectionDiffusionStrategy")
            .def( init<	ModelPart&, LinearSolverType::Pointer,	bool, int	>() )
            .def("Clear",&ResidualBasedSemiEulerianConvectionDiffusionStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::Clear)
            ;        

    class_< ResidualBasedConvectionDiffusionStrategyNonLinear< SparseSpaceType, LocalSpaceType, LinearSolverType >,
            ResidualBasedConvectionDiffusionStrategyNonLinear< SparseSpaceType, LocalSpaceType, LinearSolverType >::Pointer,
            BaseSolvingStrategyType >
            (m,"ResidualBasedConvectionDiffusionStrategyNonLinear")
            .def( init<	ModelPart&, LinearSolverType::Pointer,	bool, int, int ,double	>() )
            .def("Clear",&ResidualBasedConvectionDiffusionStrategyNonLinear< SparseSpaceType, LocalSpaceType, LinearSolverType >::Clear)
            ;
}

}  // namespace Python.

} // Namespace Kratos

