//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Angel Celigueta
//
//

// External includes

// Project includes
#include "custom_python/add_custom_processes_to_python.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "includes/kratos_parameters.h"

// Processes
#include "custom_processes/apply_laser_process.hpp"
#include "custom_processes/hypoelastic_solid_stress_tensor_calculate_process.h"
#include "custom_processes/fluid_pressure_calculate_process.h"
#include "custom_processes/merge_model_parts_process.h"



namespace Kratos
{

namespace Python
{

void  AddCustomProcessesToPython(pybind11::module& m)
{

    namespace py = pybind11;

    // Apply table values
    py::class_<ApplyLaserProcess, ApplyLaserProcess::Pointer, Process>
    (m, "ApplyLaserProcess")
    .def(py::init < ModelPart&, Parameters>())
    .def("ApplyLaser", &ApplyLaserProcess::ApplyLaser);

    py::class_<HypoelasticStressCalculateProcess, HypoelasticStressCalculateProcess::Pointer, Process >(m,"HypoelasticStressCalculateProcess")
    .def(py::init<ModelPart&, unsigned int>())
    ;
    
    py::class_<FluidPressureCalculateProcess, FluidPressureCalculateProcess::Pointer, Process >(m,"FluidPressureCalculateProcess")
    .def(py::init<ModelPart&, unsigned int>())
    ;
    
    py::class_<MergeModelPartsProcess, MergeModelPartsProcess::Pointer, Process >(m,"MergeModelPartsProcess")
    .def(py::init<> ())
    .def("MergeParts", &MergeModelPartsProcess::MergeParts)
    ;

}

}  // namespace Python.
} // Namespace Kratos

