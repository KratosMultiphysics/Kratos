//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Andreas Winterstein (a.winterstein@tum.de)
//

// System includes

#if defined(KRATOS_PYTHON)
// External includes

// Project includes
#include "includes/define_python.h"
#include "add_trilinos_strategies_to_python.h"

//Trilinos includes
#include "Epetra_FEVector.h"

// Project includes
#include "trilinos_space.h"
#include "spaces/ublas_space.h"

#include "solving_strategies/strategies/solving_strategy.h"
#include "linear_solvers/linear_solver.h"

// Strategies
#include "custom_strategies/trilinos_laplacian_meshmoving_strategy.h"
#include "custom_strategies/trilinos_structural_meshmoving_strategy.h"


namespace Kratos {
namespace Python {

void AddMeshMovingStrategies(pybind11::module& m)
{
    namespace py = pybind11;
    typedef TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector> TrilinosSparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> TrilinosLocalSpaceType;
    typedef LinearSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType > TrilinosLinearSolverType;

    typedef SolvingStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType > TrilinosBaseSolvingStrategyType;

    typedef TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector> TrilinosSparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> TrilinosLocalSpaceType;

    using TrilinosLaplacianMeshMovingStrategyType = TrilinosLaplacianMeshMovingStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType>;
    py::class_<TrilinosLaplacianMeshMovingStrategyType, typename TrilinosLaplacianMeshMovingStrategyType::Pointer, TrilinosBaseSolvingStrategyType>
    (m,"TrilinosLaplacianMeshMovingStrategy").def(py::init<Epetra_MpiComm&, ModelPart&, TrilinosLinearSolverType::Pointer, int, bool, bool, bool, int>());

    using TrilinosStructuralMeshMovingStrategyType = TrilinosStructuralMeshMovingStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType>;
    py::class_<TrilinosStructuralMeshMovingStrategyType, typename TrilinosStructuralMeshMovingStrategyType::Pointer, TrilinosBaseSolvingStrategyType>
    (m,"TrilinosStructuralMeshMovingStrategy").def(py::init<Epetra_MpiComm&, ModelPart&, TrilinosLinearSolverType::Pointer, int, bool, bool, bool, int>());
}

} // namespace Python.
} // namespace Kratos.

#endif // KRATOS_PYTHON defined
