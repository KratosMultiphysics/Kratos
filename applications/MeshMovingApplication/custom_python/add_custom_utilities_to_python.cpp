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
// Project includes
#include "custom_python/add_custom_utilities_to_python.h"

#include "custom_utilities/ball_vertex_meshmoving.h"
#include "custom_utilities/ball_vertex_meshmoving3D.h"
#include "custom_utilities/explicit_mesh_moving_utilities.h"
#include "custom_utilities/calculate_mesh_velocity_utility.h"
#include "linear_solvers/linear_solver.h"
#include "spaces/ublas_space.h"

namespace Kratos {
namespace Python {

void AddCustomUtilitiesToPython(pybind11::module& m) {
    namespace py = pybind11;

    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType> LinearSolverType;

    py::class_<BallVertexMeshMoving<2, SparseSpaceType, LinearSolverType> >(m,"BallVertexMeshMoving2D")
        .def(py::init<>())
        .def("ConstructSystem", &BallVertexMeshMoving<2, SparseSpaceType, LinearSolverType>::ConstructSystem)
        .def("BuildAndSolveSystem", &BallVertexMeshMoving<2, SparseSpaceType, LinearSolverType>::BuildAndSolveSystem)
        .def("ClearSystem", &BallVertexMeshMoving<2, SparseSpaceType, LinearSolverType>::ClearSystem);

    py::class_<BallVertexMeshMoving3D<3, SparseSpaceType, LinearSolverType>>(m,"BallVertexMeshMoving3D")
        .def(py::init<>())
        .def("ConstructSystem", &BallVertexMeshMoving3D<3, SparseSpaceType, LinearSolverType>::ConstructSystem)
        .def("BuildAndSolveSystem", &BallVertexMeshMoving3D<3, SparseSpaceType, LinearSolverType>::BuildAndSolveSystem)
        .def("ClearSystem", &BallVertexMeshMoving3D<3, SparseSpaceType, LinearSolverType>::ClearSystem);

    py::class_<ExplicitMeshMovingUtilities>(m,"ExplicitMeshMovingUtilities")
        .def(py::init<ModelPart&, ModelPart&, const double>())
        .def("ComputeExplicitMeshMovement",&ExplicitMeshMovingUtilities::ComputeExplicitMeshMovement)
        .def("FillVirtualModelPart",&ExplicitMeshMovingUtilities::FillVirtualModelPart)
        .def("ProjectVirtualValues2D",&ExplicitMeshMovingUtilities::ProjectVirtualValues<2>)
        .def("ProjectVirtualValues3D",&ExplicitMeshMovingUtilities::ProjectVirtualValues<3>)
        .def("UndoMeshMovement",&ExplicitMeshMovingUtilities::UndoMeshMovement);

    py::class_<CalculateMeshVelocityUtility>(m,"CalculateMeshVelocityUtility")
        .def(py::init<ModelPart&, Parameters>())
        .def("CalculateMeshVelocities",&CalculateMeshVelocityUtility::CalculateMeshVelocities)
        .def_static("GetMinimumBufferSize",&CalculateMeshVelocityUtility::GetMinimumBufferSize)
        ;
}

} // namespace Python.

} // Namespace Kratos
