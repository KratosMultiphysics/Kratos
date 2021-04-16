//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license:
// kratos/license.txt
//
//  Main authors:    Andreas Winterstein (a.winterstein@tum.de)
//

// System includes

// External includes
#include "includes/define_python.h"
#include "boost/numeric/ublas/vector.hpp"

// Project includes
#include "spaces/ublas_space.h"
#include "custom_python/add_custom_strategies_to_python.h"

// strategies
#include "custom_strategies/strategies/laplacian_meshmoving_strategy.h"
#include "custom_strategies/strategies/structural_meshmoving_strategy.h"
// linear solvers
#include "linear_solvers/linear_solver.h"

namespace Kratos {
namespace Python {

void AddCustomStrategiesToPython(pybind11::module& m) {
    namespace py = pybind11;

    typedef UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

    typedef LinearSolver<SparseSpaceType, LocalSpaceType> LinearSolverType;
    typedef SolvingStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType> BaseSolvingStrategyType;

    py::class_<LaplacianMeshMovingStrategy<SparseSpaceType, LocalSpaceType,LinearSolverType>,
        LaplacianMeshMovingStrategy<SparseSpaceType, LocalSpaceType,LinearSolverType>::Pointer,
        BaseSolvingStrategyType>(m,"LaplacianMeshMovingStrategy")
        .def(py::init<ModelPart &, LinearSolverType::Pointer, int, bool, bool, bool, int>())
        ;

    py::class_<StructuralMeshMovingStrategy<SparseSpaceType, LocalSpaceType,LinearSolverType>,
        StructuralMeshMovingStrategy<SparseSpaceType, LocalSpaceType,LinearSolverType>::Pointer,
        BaseSolvingStrategyType>(m,"StructuralMeshMovingStrategy")
        .def(py::init<ModelPart &, LinearSolverType::Pointer, int, bool, bool, bool, int, double>())
        ;
}

} // namespace Python.
} // Namespace Kratos
