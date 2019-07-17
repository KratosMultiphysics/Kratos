//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
// ==============================================================================
//  ChimeraApplication
//
//  License:         BSD License
//                   license: ChimeraApplication/license.txt
//
//  Authors:        Aditya Ghantasala, https://github.com/adityaghantasala
// 					Navaneeth K Narayanan
//					Rishith Ellath Meethal
// ==============================================================================
//
// System includes

// External includes

// Project includes
#include "processes/process.h"
#include "spaces/ublas_space.h"
#include "custom_python/add_custom_processes_to_python.h"

#include "custom_processes/apply_chimera_process.h"
#include "custom_processes/apply_chimera_process_monolithic.h"
#include "custom_processes/apply_chimera_process_fractional_step.h"
#include "custom_processes/rotate_region_process.h"
#include "custom_processes/calculate_signed_distance_to_2d_condition_skin_process.h"
#include "processes/calculate_signed_distance_to_3d_condition_skin_process.h"
#include "processes/calculate_distance_to_skin_process.h"

namespace Kratos
{

namespace Python
{

void AddCustomProcessesToPython(pybind11::module &m)
{

    namespace py = pybind11;

    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

    typedef ApplyChimera<2, SparseSpaceType, LocalSpaceType> BaseApplyChimera2D;
    typedef ApplyChimera<3, SparseSpaceType, LocalSpaceType> BaseApplyChimera3D;

    typedef ApplyChimeraProcessMonolithic<2, SparseSpaceType, LocalSpaceType> ApplyChimeraMonolithic2DType;
    typedef ApplyChimeraProcessMonolithic<3, SparseSpaceType, LocalSpaceType> ApplyChimeraMonolithic3DType;

    typedef ApplyChimeraProcessFractionalStep<2, SparseSpaceType, LocalSpaceType> ApplyChimeraFractionalStep2DType;
    typedef ApplyChimeraProcessFractionalStep<3, SparseSpaceType, LocalSpaceType> ApplyChimeraFractionalStep3DType;

    py::class_<BaseApplyChimera2D, BaseApplyChimera2D::Pointer, Process>(m, "BaseApplyChimera2D")
        .def(py::init<ModelPart &, Parameters>())
        .def("SetEchoLevel", &BaseApplyChimera2D::SetEchoLevel) ;
    py::class_<BaseApplyChimera3D, BaseApplyChimera3D::Pointer, Process>(m, "BaseApplyChimera3D")
        .def(py::init<ModelPart &, Parameters>())
        .def("SetEchoLevel", &BaseApplyChimera3D::SetEchoLevel) ;

    py::class_<ApplyChimeraMonolithic2DType, ApplyChimeraMonolithic2DType::Pointer, BaseApplyChimera2D>(m, "ApplyChimeraProcessMonolithic2d")
        .def(py::init<ModelPart &, Parameters>());

    py::class_<ApplyChimeraMonolithic3DType, ApplyChimeraMonolithic3DType::Pointer, BaseApplyChimera3D>(m, "ApplyChimeraProcessMonolithic3d")
        .def(py::init<ModelPart &, Parameters>());

    py::class_<ApplyChimeraFractionalStep2DType, ApplyChimeraFractionalStep2DType::Pointer, BaseApplyChimera2D>(m, "ApplyChimeraProcessFractionalStep2d")
        .def(py::init<ModelPart &, Parameters>());

    py::class_<ApplyChimeraFractionalStep3DType, ApplyChimeraFractionalStep3DType::Pointer, BaseApplyChimera3D>(m, "ApplyChimeraProcessFractionalStep3d")
        .def(py::init<ModelPart &, Parameters>());

    py::class_<RotateRegionProcess, RotateRegionProcess::Pointer, Process>(m, "RotateRegionProcess")
        .def(py::init<ModelPart &, Parameters>())
        .def("SetAngularVelocity", &RotateRegionProcess::SetAngularVelocity);
}

} // namespace Python.

} // Namespace Kratos
