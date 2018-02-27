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
#include "spaces/ublas_space.h"
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

// Project includes
#include "custom_python/add_custom_strategies_to_python.h"

// strategies
#include "custom_strategies/strategies/laplacian_meshmoving_strategy.h"
#include "custom_strategies/strategies/structural_meshmoving_strategy.h"
// linear solvers
#include "linear_solvers/linear_solver.h"

namespace Kratos {

namespace Python {
using namespace boost::python;

void AddCustomStrategiesToPython() {
  typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
  typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

  typedef LinearSolver<SparseSpaceType, LocalSpaceType> LinearSolverType;
  typedef SolvingStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>
      BaseSolvingStrategyType;
  using LaplacianMeshMovingStrategyType = LaplacianMeshMovingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >;
  using StructuralMeshMovingStrategyType = StructuralMeshMovingStrategy < SparseSpaceType, LocalSpaceType, LinearSolverType >;


  class_<LaplacianMeshMovingStrategy<SparseSpaceType, LocalSpaceType,
                                     LinearSolverType>,
         bases<BaseSolvingStrategyType>, boost::noncopyable>(
      "LaplacianMeshMovingStrategy",
      init<ModelPart &, LinearSolverType::Pointer, int, bool, bool, int>())
      .def("UpdateReferenceMesh",&LaplacianMeshMovingStrategyType::UpdateReferenceMesh)
      ;

  class_<StructuralMeshMovingStrategy<SparseSpaceType, LocalSpaceType,
                                      LinearSolverType>,
         bases<BaseSolvingStrategyType>, boost::noncopyable>(
      "StructuralMeshMovingStrategy",
      init<ModelPart &, LinearSolverType::Pointer, int, bool, bool, int>())
      .def("UpdateReferenceMesh",&StructuralMeshMovingStrategyType::UpdateReferenceMesh)
      ;
}

} // namespace Python.

} // Namespace Kratos
