//
// ==============================================================================
//  ChimeraApplication
//
//  License:         BSD License
//                   license: ChimeraApplication/license.txt
//
//  Main authors:    Aditya Ghantasala, https://github.com/adityaghantasala
//                   Navaneeth K Narayanan
//
// ==============================================================================

// System includes

// External includes

// Project includes
#include "processes/process.h"
#include "includes/model_part.h"
#include "custom_python/add_custom_processes_to_python.h"

#include "custom_processes/apply_chimera_process_monolithic.h"
#include "custom_processes/apply_chimera_process_fractional_step.h"
#include "custom_processes/rotate_region_process.h"
#include "custom_processes/calculate_signed_distance_to_2d_condition_skin_process.h"
#include "processes/calculate_signed_distance_to_3d_condition_skin_process.h"

namespace Kratos
{

namespace Python
{

void AddCustomProcessesToPython(pybind11::module &m)
{

    namespace py = pybind11;
    typedef CalculateSignedDistanceTo2DConditionSkinProcess DistanceCalculator2DType;
    typedef CalculateSignedDistanceTo3DConditionSkinProcess DistanceCalculator3DType;

    typedef ApplyChimeraProcessMonolithic<2, DistanceCalculator2DType> ApplyChimeraMonolithic2DType;
    typedef ApplyChimeraProcessMonolithic<3, DistanceCalculator3DType> ApplyChimeraMonolithic3DType;

    typedef ApplyChimeraProcessFractionalStep<2, DistanceCalculator2DType> ApplyChimeraFractionalStep2DType;
    typedef ApplyChimeraProcessFractionalStep<3, DistanceCalculator3DType> ApplyChimeraFractionalStep3DType;

    py::class_<ApplyChimeraMonolithic2DType, ApplyChimeraMonolithic2DType::Pointer, Process>(m, "ApplyChimeraProcessMonolithic2d")
        .def(py::init<ModelPart &, Parameters>());

    py::class_<ApplyChimeraMonolithic3DType, ApplyChimeraMonolithic3DType::Pointer, Process>(m, "ApplyChimeraProcessMonolithic3d")
        .def(py::init<ModelPart &, Parameters>());

    py::class_<ApplyChimeraFractionalStep2DType, ApplyChimeraFractionalStep2DType::Pointer, ApplyChimeraMonolithic2DType>(m, "ApplyChimeraProcessFractionalStep2d")
        .def(py::init<ModelPart &, Parameters>());

    py::class_<ApplyChimeraFractionalStep3DType, ApplyChimeraFractionalStep3DType::Pointer, ApplyChimeraMonolithic3DType>(m, "ApplyChimeraProcessFractionalStep3d")
        .def(py::init<ModelPart &, Parameters>());

    py::class_<RotateRegionProcess, RotateRegionProcess::Pointer, Process>(m, "RotateRegionProcess")
        .def(py::init<ModelPart &, Parameters>())
        .def("SetAngularVelocity", &RotateRegionProcess::SetAngularVelocity);
}

} // namespace Python.

} // Namespace Kratos
