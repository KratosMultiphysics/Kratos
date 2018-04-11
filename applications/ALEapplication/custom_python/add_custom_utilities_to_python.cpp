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
#include <boost/python.hpp>

// Project includes
#include "custom_python/add_custom_utilities_to_python.h"
#include "processes/process.h"

#include "custom_utilities/ball_vertex_meshmoving.h"
#include "custom_utilities/ball_vertex_meshmoving3D.h"
#include "linear_solvers/linear_solver.h"
#include "spaces/ublas_space.h"

namespace Kratos {

namespace Python {

void AddCustomUtilitiesToPython() {
  using namespace boost::python;

  typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
  typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
  typedef LinearSolver<SparseSpaceType, LocalSpaceType> LinearSolverType;

  class_<BallVertexMeshMoving<2, SparseSpaceType, LinearSolverType>,
         boost::noncopyable>("BallVertexMeshMoving2D", init<>())
      .def("ConstructSystem",
           &BallVertexMeshMoving<2, SparseSpaceType,
                                 LinearSolverType>::ConstructSystem)
      .def("BuildAndSolveSystem",
           &BallVertexMeshMoving<2, SparseSpaceType,
                                 LinearSolverType>::BuildAndSolveSystem)
      .def("ClearSystem", &BallVertexMeshMoving<2, SparseSpaceType,
                                                LinearSolverType>::ClearSystem);

  class_<BallVertexMeshMoving3D<3, SparseSpaceType, LinearSolverType>,
         boost::noncopyable>("BallVertexMeshMoving3D", init<>())
      .def("ConstructSystem",
           &BallVertexMeshMoving3D<3, SparseSpaceType,
                                   LinearSolverType>::ConstructSystem)
      .def("BuildAndSolveSystem",
           &BallVertexMeshMoving3D<3, SparseSpaceType,
                                   LinearSolverType>::BuildAndSolveSystem)
      .def("ClearSystem",
           &BallVertexMeshMoving3D<3, SparseSpaceType,
                                   LinearSolverType>::ClearSystem);
}

} // namespace Python.

} // Namespace Kratos
