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
#include "boost/numeric/ublas/vector.hpp"

// Project includes
#include "linear_solvers/linear_solver.h"
#include "spaces/ublas_space.h"

// Application includes
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/fixed_mesh_ale_utilities.h"
#include "custom_utilities/mesh_velocity_calculation.h"
#include "custom_utilities/move_mesh_utilities.h"
#include "custom_utilities/affine_transform.h"
#include "custom_utilities/parametric_affine_transform.h"

namespace Kratos {
namespace Python {

void AddCustomUtilitiesToPython(pybind11::module& m) {
    namespace py = pybind11;

    py::class_<FixedMeshALEUtilities, FixedMeshALEUtilities::Pointer>(m, "FixedMeshALEUtilities")
        .def(py::init<Model &, Parameters &>())
        .def("Initialize", &FixedMeshALEUtilities::Initialize)
        .def("SetVirtualMeshValuesFromOriginMesh", &FixedMeshALEUtilities::SetVirtualMeshValuesFromOriginMesh)
        .def("ComputeMeshMovement", &FixedMeshALEUtilities::ComputeMeshMovement)
        .def("ProjectVirtualValues2D", &FixedMeshALEUtilities::ProjectVirtualValues<2>)
        .def("ProjectVirtualValues3D", &FixedMeshALEUtilities::ProjectVirtualValues<3>)
        .def("UndoMeshMovement", &FixedMeshALEUtilities::UndoMeshMovement);

    py::class_<AffineTransform, AffineTransform::Pointer>(m, "AffineTransform")
        .def(py::init<const array_1d<double,3>&, const double, const array_1d<double,3>&, const array_1d<double,3>&>())
        .def(py::init<const array_1d<double,3>&, const array_1d<double,3>&, const array_1d<double,3>&>())
        .def(py::init<const Quaternion<double>&, const array_1d<double,3>&, const array_1d<double,3>&>())
        .def("SetRotation", static_cast<void (AffineTransform::*)(const array_1d<double,3>&, const double, const array_1d<double,3>&)>(&AffineTransform::SetRotation))
        .def("SetRotation", [](AffineTransform& rThis, const array_1d<double,3>& rAxis, const double angle, const array_1d<double,3>& rReferencePoint){rThis.SetRotation(rAxis, angle, rReferencePoint);})
        .def("SetRotation", [](AffineTransform& rThis, const array_1d<double,3>& rEulerAngles, const array_1d<double,3>& rReferencePoint){rThis.SetRotation(rEulerAngles, rReferencePoint);})
        .def("SetRotation", [](AffineTransform& rThis, const Quaternion<double>& rQuaternion, const array_1d<double,3>& rReferencePoint){rThis.SetRotation(rQuaternion, rReferencePoint);})
        .def("SetTranslation", &AffineTransform::SetTranslation)
        .def("Apply", &AffineTransform::Apply)
        ;

    py::class_<ParametricAffineTransform, ParametricAffineTransform::Pointer>(m, "ParametricAffineTransform")
        .def(py::init<const Parameters, const Parameters, const Parameters, const Parameters>())
        .def(py::init<const Parameters, const Parameters, const Parameters>())
        .def("Apply", &ParametricAffineTransform::Apply)
        ;

    void (*CalculateMeshVelocitiesBDF1)(ModelPart&, const TimeDiscretization::BDF1&) = &MeshVelocityCalculation::CalculateMeshVelocities;
    void (*CalculateMeshVelocitiesBDF2)(ModelPart&, const TimeDiscretization::BDF2&) = &MeshVelocityCalculation::CalculateMeshVelocities;
    void (*CalculateMeshVelocitiesNewmark)(ModelPart&, const TimeDiscretization::Newmark&) = &MeshVelocityCalculation::CalculateMeshVelocities;
    void (*CalculateMeshVelocitiesBossak)(ModelPart&, const TimeDiscretization::Bossak&) = &MeshVelocityCalculation::CalculateMeshVelocities;
    void (*CalculateMeshVelocitiesGeneralizedAlpha)(ModelPart&, const TimeDiscretization::GeneralizedAlpha&) = &MeshVelocityCalculation::CalculateMeshVelocities;

    m.def("CalculateMeshVelocities", CalculateMeshVelocitiesBDF1 );
    m.def("CalculateMeshVelocities", CalculateMeshVelocitiesBDF2 );
    m.def("CalculateMeshVelocities", CalculateMeshVelocitiesNewmark );
    m.def("CalculateMeshVelocities", CalculateMeshVelocitiesBossak );
    m.def("CalculateMeshVelocities", CalculateMeshVelocitiesGeneralizedAlpha );

    m.def("MoveMesh", &MoveMeshUtilities::MoveMesh );
    m.def("MoveModelPart", static_cast<void(*)(ModelPart&, const array_1d<double,3>&, const double, const array_1d<double,3>&, const array_1d<double,3>&)>(&MoveMeshUtilities::MoveModelPart));
    m.def("MoveModelPart", static_cast<void(*)(ModelPart&, const Parameters, const Parameters, const Parameters, const Parameters)>(&MoveMeshUtilities::MoveModelPart));
    m.def("MoveModelPart", static_cast<void(*)(ModelPart&, const AffineTransform&)>(&MoveMeshUtilities::MoveModelPart));
    m.def("MoveModelPart", static_cast<void(*)(ModelPart&, ParametricAffineTransform&)>(&MoveMeshUtilities::MoveModelPart));
    m.def("SuperImposeMeshDisplacement", &MoveMeshUtilities::SuperImposeMeshDisplacement );
    m.def("SuperImposeMeshVelocity", &MoveMeshUtilities::SuperImposeMeshVelocity);

}

} // namespace Python.

} // Namespace Kratos
